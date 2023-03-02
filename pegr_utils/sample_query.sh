#!/bin/bash

APIKEY=$(cat $HOME/.tokens/pegr_token.txt)
JSON=$(curl -X POST -H "Content-Type: application/json" -d "{\"id\": $1, \"userEmail\": \"jsc6015@psu.edu\"}" "https://thanos.vmhost.psu.edu/pegr/api/fetchSampleData?apiKey=$APIKEY" 2> /dev/null)
if [[ $(echo $JSON | jq -r .message) != "Success!" ]]; then
  echo "Error: $(echo $JSON | jq .message)"
  exit 1
fi

echo "ID: $(echo $JSON | jq .data[0].id)"
echo "Target: $(echo $JSON | jq .data[0].target)"
echo "Antibody: $(echo $JSON | jq .data[0].antibody)"
echo "Strain: $(echo $JSON | jq .data[0].strain)"
echo "Genetic modification: $(echo $JSON | jq .data[0].geneticModification)"
echo "Growth media: $(echo $JSON | jq .data[0].growthMedia)"
echo "Treatments: $(echo $JSON | jq .data[0].treatments)"
echo "Assay: $(echo $JSON | jq .data[0].assay)"
echo "Dedup reads: $(echo $JSON | jq .data[0].experiments[0].alignments[0].dedupUniquelyMappedReads)"
