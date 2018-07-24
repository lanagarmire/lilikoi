// npm install sax
// npm install xml-flow

const fs = require('fs')
const flow = require('xml-flow')
const { execSync } = require('child_process')

const url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
const filename = 'hmdb_metabolites.zip'
const xmlFilename = 'hmdb_metabolites.xml'
const outputFile = 'temp/hmdb_metabolites.tsv'

const downloadAndExtractFile = () => execSync(`wget ${url}; unzip ${filename}`)

const parseHmdbXmlFile = () => {
	const writeStream = fs.createWriteStream(outputFile, { flags: 'a' })
	let index = 1
	const handleMetabolite = metabolite => {
		const { accession, name, pathways } = metabolite
		let pathwaysArray = pathways ? (pathways.map ? pathways : [pathways]) : []
		const lines = pathwaysArray
			.map(pathway =>
				[accession, name, pathway.name, pathway.smpdb_id, pathway.kegg_map_id].join('\t')
			)
			.replace(/[^\x00-\x7F]/g, "") // Strip non-ascii characters b/c R sucks at UTF8
			.join('\n')

		if (lines) writeStream.write(lines + '\n')
		if (index % 1000 == 0) console.log(`${index} metabolites processed`)
		index++
	}

	xmlStream = flow(fs.createReadStream(xmlFilename))
	xmlStream.on('tag:metabolite', handleMetabolite)
	xmlStream.on('end', () => writeStream.end())
}

parseHmdbXmlFile()
