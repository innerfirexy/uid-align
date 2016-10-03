# extract subtrees from the ``parsedFull`` column in ``entropy_disf`` table from the swbd db
# Yang Xu
# 10/2/2016

import MySQLdb
import sys
import pickle
import itertools
import re

from nltk.tree import *
from nltk.probability import FreqDist


# get db connection
def db_conn(db_name):
    # db init: ssh yvx5085@brain.ist.psu.edu -i ~/.ssh/id_rsa -L 1234:localhost:3306
    conn = MySQLdb.connect(host = "127.0.0.1",
                    user = "yang",
                    port = 1234,
                    passwd = "05012014",
                    db = db_name)
    return conn

# get all sub-rules of a parsing tree
def subrules(tree_str, exclude = None):
    """
    tree_str: a str, parsing tree of a sentence
    exclude: a list of str,
    """
    if tree_str is not None:
        try:
            tree = Tree.fromstring(tree_str)
        except Exception as e:
            print('tree_str: ' + str(tree_str))
            raise e
        subtrees = list(tree.subtrees(filter = lambda t: t.height() > 2))
        rules = [str(s.productions()[0]) for s in subtrees]
        # replace all ',' and '.' with blank character
        rules = [re.sub(r'[,|\.]\s*', '', r).strip() for r in rules]
    else:
        rules = []
    # exclude
    if exclude is not None:
        rules = [r for r in rules if r not in exclude]
    return rules

# find adjacent pairs (adjacent turnIDs and different speakers)
def find_adj_pairs(data_file):
    assert isinstance(data_file, str)
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # select all convIDs
    sql = 'select distinct convID from entropy_disf'
    cur.execute(sql)
    conv_ids = [item[0] for item in cur.fetchall()]
    # within each convID, find the adjacent pairs
    result = {c_id: [] for c_id in conv_ids}
    for i, c_id in enumerate(conv_ids):
        sql = 'select distinct turnID, speaker from entropy_disf where convID = %s'
        cur.execute(sql, [c_id])
        tmp_data = cur.fetchall()
        for j in range(len(tmp_data)-1):
            if tmp_data[j][0] == tmp_data[j+1][0] - 1 and tmp_data[j][1] != tmp_data[j+1][1]:
                result[c_id].append((tmp_data[j][0], tmp_data[j+1][0]))
        # print
        sys.stdout.write('\r{}/{} convID counted'.format(i+1, len(conv_ids)))
        sys.stdout.flush()
    # save result
    pickle.dump(result, open(data_file, 'wb'))

# count repeated subrules between adjacent speakers, using the pair info saved in data_file
def count_rules(input_file, output_file):
    assert isinstance(input_file, str)
    assert isinstance(output_file, str)
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # load data file
    meta_info = pickle.load(open(input_file, 'rb'))
    # process each entry (conv_id) in meta_info
    count_dict = FreqDist()
    for i, c_id in enumerate(meta_info.keys()):
        for pair in meta_info[c_id]:
            # select all sentences for prime (pair[0])
            sql = 'select parsedFull from entropy_disf where convID = %s and turnID = %s'
            cur.execute(sql, (c_id, pair[0]))
            prime_rules = map(subrules, [item[0] for item in cur.fetchall()])
            prime_rules = list(itertools.chain.from_iterable(prime_rules))
            # select all sentences for target (pair[1])
            sql = 'select parsedFull from entropy_disf where convID = %s and turnID = %s'
            cur.execute(sql, (c_id, pair[1]))
            target_rules = map(subrules, [item[0] for item in cur.fetchall()])
            target_rules = list(itertools.chain.from_iterable(target_rules))
            # update count_dict
            for r in set(prime_rules) & set(target_rules):
                count_dict[r] += 1
        # print
        sys.stdout.write('\r{}/{} convID counted'.format(i+1, len(meta_info.keys())))
        sys.stdout.flush()
    # save to file
    pickle.dump(count_dict, open(output_file, 'wb'))

# reformat the count results so that it is readable in R
def reformat_count_result(input_file, output_file):
    assert isinstance(input_file, str)
    assert isinstance(output_file, str)
    # load
    count = pickle.load(open(input_file, 'rb'))
    # write
    with open(output_file, 'w') as fw:
        for key, val in count.items():
            fw.write(key + ', ' + str(val) + '\n')

# process all parsedFull entries and extract subrules for each
def extract_subrules_all():
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # select all data
    sql = 'select convID, globalID, parsedFull from entropy_disf'
    cur.execute(sql)
    data = cur.fetchall()
    # process each row
    for i, row in enumerate(data):
        c_id, g_id, parsed_str = row
        rules = subrules(parsed_str)
        rules_str = '~~~+~~~'.join(rules)
        # update
        sql = 'update entropy_disf set subRules = %s where convID = %s and globalID = %s'
        cur.execute(sql, (rules_str, c_id, g_id))
        # print
        sys.stdout.write('\r{}/{} row processed and updated'.format(i+1, len(data)))
        sys.stdout.flush()
    conn.commit()


# main
if __name__ == '__main__':
    # find_adj_pairs('swbd_dysf_adjacent_pairs.pkl')
    # count_rules('swbd_dysf_adjacent_pairs.pkl', 'swbd_dysf_adjacent_pairs_count.pkl')
    # reformat_count_result('swbd_dysf_adjacent_pairs_count.pkl', 'swbd_dysf_adjacent_pairs_count.txt')
    extract_subrules_all()
