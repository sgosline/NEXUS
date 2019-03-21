#!/usr/bin/env cwl-runner
class: CommandLineTool
id: salmon-quant-tool
label: salmon-quant-tool
cwlVersion: v1.0

requirements:
  - class: DockerRequirement
    dockerPull: combinelab/salmon
#  - class: InitialWorkDirRequirement
#    listing: $(inputs.output)

baseCommand: [salmon, quant, -l, A, --validateMappings]

inputs:
  mates1:
    type: File[]
    inputBinding:
      position: 1
      prefix: '-1'
  mates2:
    type: File[]
    inputBinding:
      position: 2
      prefix: '-2'
  index-dir:
    type: Directory
    inputBinding:
      prefix: -i
      position: 3
  output:
    type: string
    inputBinding:
      prefix: --output
      position: 4

outputs:
  quants:
    type: File
    outputBinding:
     glob: "*/*.sf"
  dirname:
    type: Directory
    outputBinding:
      glob: $(inputs.output)

