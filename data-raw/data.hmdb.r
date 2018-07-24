# This depends on a nodejs script.
# Run `node data.hmdb.js` first
inputFile = 'temp/hmdb_metabolites.tsv'

data.hmdb = read.delim(inputFile)
data.hmdb$smpdb_id = NULL
data.hmdb$kegg_map_id = NULL
use_data(data.hmdb, overwrite = TRUE)
