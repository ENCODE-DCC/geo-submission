import glob
import requests
import json
import ExperimentBoiler
import geoDonorMinimiser
import geoBiosampleMinimiser
import urlparse
import sys
from time import sleep


HEADERS = {'accept': 'application/json'}
GET_HEADERS = {'accept': 'application/json'}
POST_HEADERS = {'accept': 'application/json',
                'Content-Type': 'application/json'}
#SERVER = "https://test.encodedcc.org/"
SERVER = "https://www.encodeproject.org/"

def encoded_get(url, keypair=None, frame='object', return_response=False):
    url_obj = urlparse.urlsplit(url)
    new_url_list = list(url_obj)
    query = urlparse.parse_qs(url_obj.query)
    if 'format' not in query:
        new_url_list[3] += "&format=json"
    if 'frame' not in query:
        new_url_list[3] += "&frame=%s" % (frame)
    if 'limit' not in query:
        new_url_list[3] += "&limit=all"
    if new_url_list[3].startswith('&'):
        new_url_list[3] = new_url_list[3].replace('&', '', 1)
    get_url = urlparse.urlunsplit(new_url_list)
    max_retries = 10
    max_sleep = 10
    while max_retries:
        try:
            if keypair:
                response = requests.get(get_url,
                                        auth=keypair,
                                        headers=GET_HEADERS)
            else:
                response = requests.get(get_url, headers=GET_HEADERS)
        except (requests.exceptions.ConnectionError,
                requests.exceptions.SSLError) as e:
            print >> sys.stderr, e
            sleep(max_sleep - max_retries)
            max_retries -= 1
            continue
        else:
            if return_response:
                return response
            else:
                return response.json()


def getKeyPair(path_to_key_pair_file, server_name):
    keysf = open(path_to_key_pair_file, 'r')
    keys_json_string = keysf.read()
    keysf.close()
    keys = json.loads(keys_json_string)
    key_dict = keys[server_name]
    AUTHID = key_dict['key']
    AUTHPW = key_dict['secret']
    return (AUTHID, AUTHPW)


def extract_biosamples(exp):
    samples = []
    if exp['status'] == 'released' and \
       'replicates' in exp and \
       len(exp['replicates']) > 0:
        for replicate in exp['replicates']:
            if replicate['status'] == 'released' and \
               replicate['library']['status'] == 'released' and \
               replicate['library']['biosample']['status'] == 'released':
                samples.append(replicate['library']['biosample']['accession'])
    return list(set(samples))


def extract_controls(exp):
    if "possible_controls" in exp and \
       len(exp['possible_controls']) > 0:
        controls_list = []
        for e in exp['possible_controls']:
            controls_list.append(e['accession'])

        return list(set(controls_list))
    else:
        return []


def extract_donors(biosamples_list):
    donors = []
    for biosample in biosamples_list:
        if biosample['status'] == 'released' and \
           'donor' in biosample and \
           biosample['donor']['status'] == 'released':
            donors.append(biosample['donor']['accession'])
    return list(set(donors))


keypair = getKeyPair('keypairs.json', 'test')

AUTHID = keypair[0]
AUTHPW = keypair[1]

# phase 1 - collect all experiments submitted so far.

submittedExperiments = set()
for filename in glob.glob('../experiments/*.json'):
    submittedExperiments.add(filename.split('/')[2].split('_')[0])

e3 =0
other =0
m = 0
f_e3 = open('e3_submitted_to_geo.tsv', "w")
x = open('not_e3_submitted_to_geo.tsv', "w")

for experiment in submittedExperiments:
    URL = SERVER + experiment + "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    experiment_o = response.json()
    if experiment_o['award']['rfa']=='ENCODE3':
        e3 += 1
        f_e3.write(experiment + "\t" + str(experiment_o['dbxrefs']) + '\t' +experiment_o['award']['rfa'] + '\n')
    else:
        other += 1
        x.write(experiment + "\t" + str(experiment_o['dbxrefs']) + '\t' + experiment_o['award']['rfa']+ '\n')
    m += 1
    if m % 10 == 0:
        print ('processed ' + str(m))

print ('E3 = ' + str(e3) + '  other = ' + str(other))
f_e3.close()
x.close()
