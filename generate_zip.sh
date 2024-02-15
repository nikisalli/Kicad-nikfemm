#!/bin/bash

rm NikFEMM.zip
mv metadata.json metadata_.json
jq --arg today "$(date +%Y-%m-%d)" '.versions[0].version |= $today' metadata_.json > metadata.json

git ls-files  -- 'python-kicad-plugin/metadata.json' 'python-kicad-plugin/resources*.png' 'python-kicad-plugin/plugins*.png' 'python-kicad-plugin/plugins*.py' | xargs zip NikFEMM.zip
mv metadata_.json metadata.json