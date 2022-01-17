import requests
import click
import csv
import json 

rq_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"

tmplt = {"PRIDE dataset URL": 'N/A',
            "Lab Head": 'N/A',
            "Correspondence": 'N/A',  # have mercy
            "Affiliation": 'N/A',
            "Original dataset submitter": 'N/A',
            "E-mail": 'N/A',
            "PubMed ID": 'N/A',
            "Project description": 'N/A',
            "Project title": 'N/A'}

@click.command()
@click.option('--accession', '-a', type=str,
                help="PRIDE database accessions start with PXD. You must make sure yourself that it exists, otherwise this will fail.")
@click.option('--outputfilename', '-o', type=click.Path(writable=True), default="REPLACEME",
                help='Will default to the given accession or NA.')
def query_PRIDEAPI(accession: str, outputfilename: click.Path) -> None:
    response = requests.get(rq_url.format(accession=accession))
    if response.ok:
        #print(json.dumps(response.json(), indent=4, sort_keys=True))
        tmplt['PubMed ID'] = response.json()["references"][0]["pubmedId"]
        tmplt['PRIDE dataset URL'] = rq_url.format(accession=accession)

        tmplt['Lab Head'] = ' '.join([response.json()["labPIs"][0]["title"],response.json()["labPIs"][0]["name"]])
        tmplt['Affiliation'] = response.json()["labPIs"][0]["affiliation"]
        tmplt['Correspondence'] = response.json()["labPIs"][0]["email"]
        
        tmplt['Project description'] = response.json()["projectDescription"]
        tmplt['Project title'] = response.json()["title"]

        tmplt['Original dataset submitter'] = ' '.join([response.json()["submitters"][0]["title"],response.json()["submitters"][0]["name"]])
        tmplt['E-mail'] = response.json()["submitters"][0]["email"]
        
    else:
        Warning("Could not retrieve PXD meta information from PRIDE-API.")

    if outputfilename == "REPLACEME":
        if not accession:
                with click.get_current_context() as ctx:
                        click.echo(ctx.get_help())
                raise click.Abort

        gen = ''.join([accession, '.meta'])
        #click.Path(gen, writable=True)
        fin(gen)
    else:     
        fin(str(outputfilename))

def fin(on: str) -> None:
    #'Will default to the given accession or NA.'
    with open(on, 'w') as f:
        w = csv.DictWriter(f, tmplt.keys())
        w.writeheader()
        w.writerow(tmplt)

if __name__ == '__main__':
    query_PRIDEAPI()
