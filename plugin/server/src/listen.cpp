
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <string.h>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/monte_carlo.h"
#include "feasst/include/feasst.h"
#include "server/include/listen.h"

using namespace std;

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
  int code;
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
    code = read(client_socket, buffer, buffer_size_/2);
    ASSERT(code >= 0, "error");
    if (strcmp(buffer, "EndListen") == 0) {
      INFO(buffer);
      finished = true;
    } else {
      strcpy(buffer, "test");
      code = write(client_socket, buffer, strlen(buffer));
    }
    ASSERT(code >= 0, "error");
  }
  close(server_socket);
  delete[] buffer;
}

}  // namespace feasst
