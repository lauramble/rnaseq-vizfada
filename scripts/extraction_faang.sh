#!/usr/bin/sh

echo $1

# extraction de toutes les données FAANG Sus scrofa
curl -POST "http://data.faang.org/api/file/_search/?size=25000" -d '{
    "query": {
      "bool": {
        "must": [
          {"match": {"species.text": "'"$1"'"}},
          {"match": {"experiment.standardMet": "FAANG"}}
        ]
      }
    }  
}' > species.json

# extraction de toutes les données FAANG specimen
curl -POST "http://data.faang.org/api/specimen/_search/?size=20000" -d '{
  "query": {
    "bool": {
        "must":
          {"match": {"standardMet": "FAANG"}}
    }
  }
}' > specimens.json

# extraction de toutes les données FAANG RNA-seq experiment
curl -POST "http://data.faang.org/api/experiment/_search/?size=20000" -d '{
    "_source": "_id",
    "query": {
        "bool": {
            "must": [{
                "exists": {"field": "RNA-seq"}                
                },{
                "match" : {"standardMet": "FAANG"}
                }]
            
        }
    }  
}' > experiments.json

Rscript GetMeta.R specimens.json experiments.json species.json
