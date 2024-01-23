
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <string.h>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/monte_carlo.h"
#include "feasst/include/feasst.h"
#include "server/include/listen.h"

//using namespace std;

namespace feasst {

Listen::Listen(argtype * args) {
  class_name_ = "Listen";
  port_ = integer("port", args, 50007);
  buffer_size_ = integer("buffer_size", args, 1000);
  ASSERT(buffer_size_ % 2 == 0, "buffer_size must be even.");
}
Listen::Listen(argtype args) : Listen(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapListen {
 public:
  MapListen() {
    auto obj = MakeListen();
    obj->deserialize_map()["Listen"] = obj;
  }
};

static MapListen mapper_Listen = MapListen();

Listen::Listen(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3204 && version <= 3204, "mismatch version: " << version);
  feasst_deserialize(&port_, istr);
  feasst_deserialize(&buffer_size_, istr);
}

void Listen::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3204, ostr);
  feasst_serialize(port_, ostr);
  feasst_serialize(buffer_size_, ostr);
}

// implementation from https://stackoverflow.com/questions/20732980/how-to-use-socket-with-a-python-client-and-a-c-server
void Listen::run(MonteCarlo * mc) {
  char* buffer = new char[buffer_size_ + 1];
  int size;
  int server_socket=socket(AF_INET, SOCK_STREAM, 0);
  sockaddr_in serverAddr;
  serverAddr.sin_family = AF_INET;
  serverAddr.sin_port = htons(port_);
  serverAddr.sin_addr.s_addr = INADDR_ANY;
  bind(server_socket, (struct sockaddr*)&serverAddr, sizeof(struct sockaddr));
  listen(server_socket, 1);
  sockaddr_in clientAddr;
  socklen_t sin_size = sizeof(struct sockaddr_in);
  int client_socket = accept(server_socket, (struct sockaddr*)&clientAddr, &sin_size);
  bool finished = false;
  while (!finished) {
    bzero(buffer, buffer_size_);
    size = read(client_socket, buffer, buffer_size_/2);
    DEBUG("received: " << buffer);
    DEBUG("receive size: " << size);
    ASSERT(size >= 0, "error");
    if (strcmp(buffer, "EndListen") == 0) {
      DEBUG(buffer);
      finished = true;
    } else {
      if (size > 0) {
        std::string line(buffer);
        std::pair<std::string, argtype> marg = parse_line(line, NULL, NULL);
        DEBUG(marg.first);
        DEBUG(str(marg.second));
        arglist list(1);
        list[0] = marg;
        mc->parse_args(&list);
        ASSERT(list.size() == 0, "Unrecognized argument: " << list.begin()->first);
      }
      strcpy(buffer, "Message received. Continue.");
      size = write(client_socket, buffer, strlen(buffer));
      DEBUG("sent: " << buffer);
      DEBUG("send size: " << size);
    }
    ASSERT(size >= 0, "error");
  }
  close(server_socket);
  delete[] buffer;
}

}  // namespace feasst
