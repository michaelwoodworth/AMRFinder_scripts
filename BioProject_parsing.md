# BioProject Parsing
Finding SRA accession IDs from other BioProject metadata.

## Motivation
A collaborator or publication may list a BioProject metadata identifier that is not an SRA accession, which can make command line downloads more difficult.  This code is a workaround for pulling SRA accessions with a BioProject ID.

## Workflow
1. Test esearch with BioSample query

```console
esearch -db biosample -q ${ID} | efetch
```

2. Pipe to pull SRA accession

```console
esearch -db biosample -q ${ID} | efetch | grep "SRA" | cut -d":" -f5
```

3. Loop with list of IDs and concatenate to list

```console
# loop to make list of accessions from IDs
for ID in `cat IDlist.txt`; do SRA=$(esearch -db biosample -q $ID | efetch | grep "SRA" | cut -d":" -f5); echo ${ID}, $SRA >> SRAlist.csv; echo $SRA added; done

# loop to second list of values without missing second field (SRA accession)
for test in `cat SRAlist.csv | cut -d, -f2`; do echo $test >> SRA2.txt; done
```

4. Save file with IDs & successful accession number searches

```console
for accession in `cat SRA2.txt`; do grep $accession SRAlist.csv >> SRA3.csv; done
```