#!python3
# extract subtrees from the ``parsedFull`` column in ``entropy_disf`` table from the swbd db
# Yang Xu
# 10/2/2016

import MySQLdb
import sys
import pickle
import itertools
import re
import math

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

# obtain the frequency of subrules
def subrules_freq(output_file):
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # select all data
    sql = 'select subRules from entropy_disf where subRules <> \'\'' # non-empty
    cur.execute(sql)
    data = cur.fetchall()
    # count freq
    count = FreqDist()
    for item in data:
        rules = item[0].split('~~~+~~~')
        for r in rules:
            count[r] += 1
    # save to R-friendly format
    with open(output_file, 'w') as fw:
        for key, val in count.items():
            fw.write(key + ', ' + str(val) + '\n')

# create a dict that has (convID, turnID) as keys, and all subrules within as values
# all turnIDs from the entropy_disf table are included
def turn_subrules_dict(output_file):
    assert isinstance(output_file, str)
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # select all unique (convID, turnID) pairs from table
    sql = 'select distinct convID, turnID from entropy_disf' # where subRules <> \'\'
    cur.execute(sql)
    all_keys = [(row[0], row[1]) for row in cur.fetchall()]
    # create dict
    result_dict = {}
    for i, key in enumerate(all_keys):
        sql = 'select subRules from entropy_disf where convID = %s and turnID = %s'
        cur.execute(sql, (key[0], key[1]))
        sents = [row[0] for row in cur.fetchall()]
        rules = []
        for s in sents:
            rules += s.split('~~~+~~~')
        result_dict[key] = rules
        # print
        sys.stdout.write('\r{}/{} keys inserted'.format(i+1, len(all_keys)))
        sys.stdout.flush()
    # save to file
    pickle.dump(result_dict, open(output_file, 'wb'))

# the function that computes the prior probability of a rule appearing in a turn
def prior(r_dict, rule):
    assert isinstance(r_dict, dict)
    assert isinstance(rule, str)
    count = 0
    for key, val in r_dict.items():
        if rule in val:
            count += 1
    return float(count) / len(r_dict)

# the function that computes the probability boost for given list of subrules
def probBoost(p_dict, r_dict, rule, method='log-odds'):
    assert isinstance(p_dict, dict)
    assert isinstance(r_dict, dict)
    assert isinstance(rule, str)
    assert method in ['log-odds', 'diff']
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # considering all (prime, target) pairs, count the number of pairs whose prime contains rule
    # and the number of pairs whose prime and target both contains rule
    count_pair = 0
    count_t = 0
    count_p = 0
    count_pt = 0
    for key, val in p_dict.items():
        for p_turn, t_turn in val:
            count_pair += 1
            if rule in r_dict[(key, t_turn)]:
                count_t += 1
            if rule in r_dict[(key, p_turn)]:
                count_p += 1
                if rule in r_dict[(key, t_turn)]:
                    count_pt += 1
    if count_p == 0:
        return None
    prob = float(count_pt) / count_p
    # compute prior probability
    # prior_prob = prior(r_dict, rule)
    prior_prob = float(count_t) / count_pair
    if prior_prob == 0:
        return None
    # calculate probability boost
    if method == 'log-odds':
        if prior_prob != 0:
            return math.log(prob / prior_prob)
        else:
            return float('inf')
    else:
        return prob - prior_prob

# count the number of subrules in specified pairs


# experiment that computes the probability boost
def exp_probBoost():
    # load pairs_dict and rules_dict
    pairs_dict = pickle.load(open('swbd_dysf_adjacent_pairs.pkl', 'rb'))
    rules_dict = pickle.load(open('turn_subrules_dict.pkl', 'rb'))
    # target rules
    target_rules = ['NP -> PRP', 'S -> NP VP', 'PP -> IN NP', 'ADVP -> RB', 'S -> VP', 'SBAR -> S', 'NP -> DT NN', 'NP -> NN', 'NP -> NP PP', 'SBAR -> IN S', \
        'NP -> NP SBAR', 'VP -> TO VP', 'NP -> DT', 'NP -> NNS', 'SBAR -> WHNP S', 'VP -> MD VP', 'VP -> VB NP', 'VP -> VBP SBAR', 'VP -> VBP VP', 'S -> S CC S']
    # compute the probability boost for each rule
    pb_results = []
    for i, rule in enumerate(target_rules):
        pb_results.append((rule, probBoost(pairs_dict, rules_dict, rule, method='log-odds')))
        sys.stdout.write('\r{}/{} rules computed'.format(i+1, len(target_rules)))
        sys.stdout.flush()
    # save results
    with open('pb_results_20.txt', 'w') as fw:
        for row in pb_results:
            fw.write(row[0] + ', ' + str(row[1]) + '\n')

# experiment that counts the number of pairs, in which a specified rule only occurs once in prime and target
def exp_unq_occur():
    pass



# main
if __name__ == '__main__':
    # find_adj_pairs('swbd_dysf_adjacent_pairs.pkl')
    # count_rules('swbd_dysf_adjacent_pairs.pkl', 'swbd_dysf_adjacent_pairs_count.pkl')
    # reformat_count_result('swbd_dysf_adjacent_pairs_count.pkl', 'swbd_dysf_adjacent_pairs_count.txt')
    # extract_subrules_all()
    # subrules_freq('all_subrules_freq.txt')
    # turn_subrules_dict('turn_subrules_dict.pkl')
    exp_probBoost()
