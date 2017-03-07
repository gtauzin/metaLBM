from __future__ import print_function

import numpy as np
from string import Template
import sys
import json
import os

def GetInputOrDefault(d, key, default):
   if key in d.keys():
      return d[key]
   else:
      return default

def GetkPairs(kmin, kmax):
   kX = []
   kY = []
   #kX = ['5']
   #kY = ['0']

   for i in range(-kmax, kmax+1):
      for j in range(-kmax, kmax+1):
         norm = i*i + j*j
         if norm >= kmin*kmin and norm <= kmax*kmax:
            kX.append(str(i))
            kY.append(str(j))
            #pass
   print(kX, kY)

   return kX, kY

def main():
   JSONfilename = sys.argv[1]

   input_folder = os.path.dirname(os.path.abspath(JSONfilename))
   print("JSONfilename: ", JSONfilename)
   with open(JSONfilename) as json_file:
      input_json = json.load(json_file)

   d_template = {}

   d_template['lengthX_g'] = input_json['lengthX']
   d_template['lengthY_g'] = input_json['lengthY']

   d_template['initDensityType'] = input_json['initDensity']
   d_template['initDensityType'] = input_json['initDensity']['type']
   d_template['initDensityValue'] = input_json['initDensity']['value']
   d_template['initVelocityType'] = input_json['initVelocity']['type']
   d_template['initVelocityXValue'] = input_json['initVelocity']['valueX']
   d_template['initVelocityYValue'] = input_json['initVelocity']['valueY']

   d_template['solverMethod'] = input_json['solverMethod']
   d_template['tau'] = input_json['tau']
   d_template['startIteration'] = input_json['startIteration']
   d_template['iterationMax'] = input_json['iterationMax']
   d_template['writeStep'] = input_json['writeStep']
   d_template['backupStep'] = input_json['backupStep']
   d_template['prefix'] = input_json['prefix']

   d_template['forcingMethod'] = input_json['forcingMethod']
   d_template['numberForces'] = len(input_json['forces'])
   forceTypeArray = []
   forceAmplitudeXArray = []
   forceAmplitudeYArray = []
   forceWaveLengthXArray = []
   forceWaveLengthYArray = []
   forceNumberWavenumberPairsArray = []
   forcekXArray = []
   forcekYArray = []
   lengthScale = 1

   forcePeriodTimeArray = []
   for k in range(0, d_template['numberForces']):
      forceTypeArray.append('ForceType::'+input_json['forces'][k]['type'])
      forceAmplitudeXArray.append(str(input_json['forces'][k]['amplitudeX']))
      forceAmplitudeYArray.append(str(input_json['forces'][k]['amplitudeY']))
      forceWaveLengthXArray.append(str(GetInputOrDefault(input_json['forces'][k], 'waveLengthX', 64)))
      forceWaveLengthYArray.append(str(GetInputOrDefault(input_json['forces'][k], 'waveLengthY', 64)))
      forcePeriodTimeArray.append(str(GetInputOrDefault(input_json['forces'][k], 'periodTime', 0.0)))
      kmin = GetInputOrDefault(input_json['forces'][k], 'kMin', 0)
      kmax = GetInputOrDefault(input_json['forces'][k], 'kMax', 0)
      tempkXArray, tempkYArray = GetkPairs(kmin, kmax)
      forceNumberWavenumberPairsArray.append(len(tempkXArray))
      forcekXArray.append('{{'+','.join(tempkXArray)+'}}')
      forcekYArray.append('{{'+','.join(tempkYArray)+'}}')
      #BEWARE THIS DEFINITION IS HIGHLY DEPENDENT ON HOW YOU FORCE
      if lengthScale < kmin + (kmax-kmin)/2:
         lengthScale = kmin + (kmax-kmin)/2

   d_template['forceTypeArray'] = ', '.join(forceTypeArray)
   d_template['forceAmplitudeXArray'] = ', '.join(forceAmplitudeXArray)
   d_template['forceAmplitudeYArray'] = ', '.join(forceAmplitudeYArray)
   d_template['forceWaveLengthXArray'] = ', '.join(forceWaveLengthXArray)
   d_template['forceWaveLengthYArray'] = ', '.join(forceWaveLengthYArray)
   d_template['forcePeriodTimeArray'] = ', '.join(forcePeriodTimeArray)
   d_template['forceNumberWavenumberPairs'] = max([0]+forceNumberWavenumberPairsArray)
   # d_template['forcekXArray'] = ', '.join(forcekXArray)
   # d_template['forcekYArray'] = ', '.join(forcekYArray)
   d_template['forcekXArray'] = '{'+', '.join(forcekXArray)+'}'
   d_template['forcekYArray'] = '{'+', '.join(forcekYArray)+'}'
   d_template['lengthScale'] = d_template['lengthX_g']/lengthScale

   d_template['numberBCs'] = len(input_json['boundaries'])
   boundaryTypeArray = []
   boundaryPositionArray = []
   boundaryStartArray = []
   boundaryEndArray = []
   boundaryPressureArray = []
   boundaryVelocityXArray = []
   boundaryVelocityYArray = []
   for k in range(0, d_template['numberBCs']):
      boundaryTypeArray.append('BoundaryType::'+input_json['boundaries'][k]['type'])
      boundaryPositionArray.append('BoundaryPosition::'+input_json['boundaries'][k]['position'])
      boundaryStartArray.append(str(input_json['boundaries'][k]['start']))
      boundaryEndArray.append(str(input_json['boundaries'][k]['end']))
      boundaryPressureArray.append(str(input_json['boundaries'][k]['pressure']))
      boundaryVelocityXArray.append(str(input_json['boundaries'][k]['velocityX']))
      boundaryVelocityYArray.append(str(input_json['boundaries'][k]['velocityY']))
   d_template['boundaryTypeArray'] = ', '.join(boundaryTypeArray)
   d_template['boundaryPositionArray'] = ', '.join(boundaryPositionArray)
   d_template['boundaryStartArray'] = ', '.join(boundaryStartArray)
   d_template['boundaryEndArray'] = ', '.join(boundaryEndArray)
   d_template['boundaryPressureArray'] = ', '.join(boundaryPressureArray)
   d_template['boundaryVelocityXArray'] = ', '.join(boundaryVelocityXArray)
   d_template['boundaryVelocityYArray'] = ', '.join(boundaryVelocityYArray)

   d_template['numberOutputs'] = len(input_json['outputs'])
   outputTypeArray = []
   for k in range(0, d_template['numberOutputs']):
      outputTypeArray.append('OutputType::'+input_json['outputs'][k]['type'])
   d_template['outputTypeArray'] = ', '.join(outputTypeArray)

   d_template['writeNextDensity'] = GetInputOrDefault(input_json, 'writeNextDensity', 1)
   d_template['writePreviousDensity'] = GetInputOrDefault(input_json, 'writePreviousDensity', 1)
   d_template['writePrevious2Density'] = GetInputOrDefault(input_json, 'writePrevious2Density', 1)

   d_template['writeNextVelocity'] = GetInputOrDefault(input_json, 'writeNextVelocity', 1)
   d_template['writePreviousVelocity'] = GetInputOrDefault(input_json, 'writePreviousVelocity', 1)
   d_template['writePrevious2Velocity'] = GetInputOrDefault(input_json, 'writePrevious2Velocity', 1)
   d_template['writeNextVorticity'] = GetInputOrDefault(input_json, 'writeNextVorticity', 1)
   d_template['writePreviousVorticity'] = GetInputOrDefault(input_json, 'writePreviousVorticity', 1)
   d_template['writeNextForce'] = GetInputOrDefault(input_json, 'writeNextForce', 1)
   d_template['writeNextAlpha'] = GetInputOrDefault(input_json, 'writeNextAlpha', 1)

   d_template['writeNextEntropy'] = GetInputOrDefault(input_json, 'writeNextEntropy', 1)
   d_template['writePreviousEntropy'] = GetInputOrDefault(input_json, 'writePreviousEntropy', 1)

   with open(input_folder +'/../inputPy/input.txt', 'r') as template_file:
      input_template = template_file.read()

   template = Template(input_template)
   input = template.safe_substitute(d_template)

   Hfilename = input_folder +'/../inputPy/input.h'

   with open(Hfilename, 'w+') as h_file:
      h_file.write(input)

if  __name__ =='__main__':
   main()
