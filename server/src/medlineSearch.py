from Bio import Entrez
from sys import stderr
from os import remove,listdir

def search_articles(query):
	Entrez.email = 'oagarwal@seas.upenn.edu'
	handle = Entrez.esearch(db='pubmed',sort='relevance',retmax='10',retmode='xml',term=query)
	id_list = Entrez.read(handle)['IdList']
  	results = fetch_details(id_list)
	return results

def fetch_details(id_list):
	ids = ','.join(id_list)
	Entrez.email = 'oagarwal@seas.upenn.edu'
	handle = Entrez.efetch(db='pubmed',retmode='xml',id=ids)
	results = Entrez.read(handle)
	article_names = [] 
	article_ids = []
	article_years = []
	for i, result in enumerate(results['PubmedArticle']):
		article_names.append(result['MedlineCitation']['Article']['ArticleTitle'])
		article_ids.append(str(result['MedlineCitation']['PMID']))
		article_years.append(result['MedlineCitation']['DateCreated']['Year'])
		save_article(str(result['MedlineCitation']['PMID']),result['MedlineCitation']['Article']['Abstract']['AbstractText'])
	results = {}
	results['names'] = article_names
	results['pmids'] = article_ids
	results['years'] = article_years
	return results

def save_article(pmid,abstract):
	abstract = '\n'.join(abstract)
	f = open("data/medline/"+pmid+".txt",'w')
	f.write(abstract.encode('utf-8'))
	f.close()

def clear_articles():
	for file_name in get_all_files("data/medline"):
		if file_name[-4:] != "conf":
			remove("data/medline/"+file_name)

def get_all_files(input_dir):
        return [direc for direc in listdir(input_dir)]

