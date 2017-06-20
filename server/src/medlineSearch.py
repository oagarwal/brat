from Bio import Entrez
from sys import stderr
from os import rename,listdir
from os.path import isdir
from config import DATA_DIR

def search_articles(query,username):
	Entrez.email = 'oagarwal@seas.upenn.edu'
	try:
		handle = Entrez.esearch(db='pubmed',sort='relevance',retmax='15',retmode='xml',term=query)
	except:
		return {}
	id_list = Entrez.read(handle)['IdList']
	if len(id_list) <=0:
		return {}
	query += " AND Randomized controlled Trial[Publication Type]"
	handle = Entrez.esearch(db='pubmed',sort='relevance',retmax='5',retmode='xml',term=query)
	id_list += Entrez.read(handle)['IdList']
  	results = fetch_details(set(id_list),username)
	return results

def fetch_details(id_list,username):
	ids = ','.join(id_list)
	Entrez.email = 'oagarwal@seas.upenn.edu'
	handle = Entrez.efetch(db='pubmed',retmode='xml',id=ids)
	results = Entrez.read(handle)
	article_names = [] 
	article_ids = []
	article_years = []
	articles_rct = []
	articles_mesh = []
	for i, result in enumerate(results['PubmedArticle']):

		if 'PMID' in result['MedlineCitation']:
			article_ids.append(str(result['MedlineCitation']['PMID']))
		else:
			continue;

                if 'MeshHeadingList' in result['MedlineCitation']:
                        mesh_terms =  set([str(val['DescriptorName']) for val in result['MedlineCitation']['MeshHeadingList']])
                        articles_mesh.append(", ".join(mesh_terms));


		if 'Article' in result['MedlineCitation']:
			if 'ArticleTitle' in result['MedlineCitation']['Article']:
				article_names.append(result['MedlineCitation']['Article']['ArticleTitle'])
			else:
				article_names.append(None)
			if 'PublicationTypeList' in result['MedlineCitation']['Article']:
				articles_rct.append('Randomized Controlled Trial' in result['MedlineCitation']['Article']['PublicationTypeList'] and 'Humans' in mesh_terms);
			else:
				articles_rct.append(None)

		if 'DateCreated' in result['MedlineCitation'] and 'Year' in result['MedlineCitation']['DateCreated']:
			article_years.append(result['MedlineCitation']['DateCreated']['Year'])
		else:
			article_years.append(None)
		
		if 'Article' in result['MedlineCitation'] and 'Abstract' in result['MedlineCitation']['Article'] and 'AbstractText' in result['MedlineCitation']['Article']['Abstract']:	
			save_article(str(result['MedlineCitation']['PMID']),result['MedlineCitation']['Article']['Abstract']['AbstractText'],username)
		else:
			save_article(str(result['MedlineCitation']['PMID']),'Abstract Unavailable',username)

	results = {}
	results['names'] = article_names
	results['pmids'] = article_ids
	results['years'] = article_years
	results['rct'] = articles_rct
	results['mesh'] = articles_mesh
	return results

def save_article(pmid,abstract,username):
	abstract = '\n'.join(abstract)
	f = open(DATA_DIR+'/'+username+'/'+pmid+".txt",'w')
	f.write(abstract.encode('utf-8'))
	f.close()

def clear_articles(username):
	for file_name in get_all_files(DATA_DIR+'/'+username):
		if file_name[-4:] != "conf":
			rename(DATA_DIR+'/'+username+'/'+file_name,DATA_DIR+'/'+username+'/old/'+file_name)
	return {}

def get_all_files(input_dir):
        return [direc for direc in listdir(input_dir) if not isdir(input_dir+'/'+direc)]

