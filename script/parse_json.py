from __future__ import print_function

import sys
import json

def main():
   JSONfilename = sys.argv[1]
   key = sys.argv[2]

   with open(JSONfilename) as json_file:
      input_json = json.load(json_file)

   print(input_json[key])

if  __name__ =='__main__':
   main()
