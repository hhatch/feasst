import json

def process(in_file_name, out_file_name):
  with open(in_file_name + ".xyz", "r") as f:
    content = f.readlines()

  out_file = open(out_file_name + ".xyz", "w")
  print(content[0], file=out_file, end='')
  print("-1 " + content[1], file=out_file, end='')

  data = {}
  atomTypes = list()
  #numTypes = list()
  for line in range(len(content)):
    if line >= 2:
      name = content[line].split(" ")[0]
      # print("name", name)
      if name not in atomTypes:
        atomTypes.append(name)
        data[name] = 1
      else:
        data[name] = data[name] + 1

  for atom in atomTypes:
    for line in range(len(content)):
      if line >= 2:
        name = content[line].split(" ")[0]
        if name == atom:
          print(content[line], file=out_file, end='')

  with open(out_file_name + "_types.json", 'w') as json_file:
    json.dump(data, json_file)
  return data

#in_file_name = "MFI_replicate"
#out_file_name = in_file_name + "_out"

data = process(in_file_name = "MFI_replicate", out_file_name = "MFI_replicate_out")

#print(content)
#print(atomTypes)
#print(data)
#print(nlines)
