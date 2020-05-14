import json

''' concatenate two xyz files assuming FEASST file format '''

in_file_name1 = "co2"
in_file_name2 = "MFI_replicate_out"
out_file_name = "both"

def read_box_bounds(content):
  lx = float(content[1].split(" ")[1])
  ly = float(content[1].split(" ")[2])
  lz = float(content[1].split(" ")[3])
  return lx, ly, lz

def print_atoms(content, out_file):
  for line in range(len(content)):
    if line >= 2:
      print(content[line], file=out_file, end='')

with open(out_file_name + ".xyz", "w") as out_file:
  with open(in_file_name1+ ".xyz", "r") as f1:
    content1 = f1.readlines()
    lx1, ly1, lz1 = read_box_bounds(content1)
    n1 = int(content1[0])
    with open(in_file_name2+ ".xyz", "r") as f2:
      content2 = f2.readlines()
      lx2, ly2, lz2 = read_box_bounds(content2)
      assert(lx1 == lx2)
      assert(ly1 == ly2)
      assert(lz1 == lz2)
      n2 = int(content2[0])

      print(str(n1 + n2) + '\n', file=out_file, end='')
      print(content1[1], file=out_file, end='')
      print_atoms(content1, out_file)
      print_atoms(content2, out_file)
