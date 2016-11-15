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
import csv

from nltk.tree import *
from nltk.probability import FreqDist
import nltk_tgrep


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
def find_adj_pairs(output_file):
    assert isinstance(output_file, str)
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
    pickle.dump(result, open(output_file, 'wb'))

# count repeated subrules between adjacent speakers, using the pair info saved in input_file
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

##
# construct the dict that maps from (conv_id, turn_id) to the list of tuples [(subrules, parsetree), ...]
# each tuple in the list corresponds to a complete sentence
# save the dict to the output_file
def subrules_parsetree_map(out_file):
    assert isinstance(out_file, str)
    # db conn
    conn = db_conn('swbd')
    cur = conn.cursor()
    # select all unique (convID, turnID) pairs from table
    sql = 'select distinct convID, turnID from entropy_disf' # where subRules <> \'\'
    cur.execute(sql)
    all_keys = [(row[0], row[1]) for row in cur.fetchall()]
    # construct dict
    result = {}
    for i, key in enumerate(all_keys):
        sql = 'select parsedFull, subRules from entropy_disf where convID = %s and turnID = %s'
        cur.execute(sql, (key[0], key[1]))
        tmp_data = cur.fetchall()
        tmp_list = []
        for row in tmp_data:
            tree_str, rules_str = row
            rules = rules_str.split('~~~+~~~')
            tmp_list.append((rules, tree_str))
        result[key] = tmp_list
        # print
        sys.stdout.write('\r{}/{} keys inserted'.format(i+1, len(all_keys)))
        sys.stdout.flush()
    # save to file
    pickle.dump(result, open(out_file, 'wb'))


##
# find the parsetree of a given subrule, and its conv_id, turn_id
# locate the list [(subrules, parsetree)] by (conv_id, turn_id), and scan through the list
# and return all the parsetrees that contains the given subrule
def find_parsetree(map_dict, conv_id, turn_id, rule):
    """
    map_dict: the resulting dict returned by subrules_parsetree_map
    return: the parse tree that contains the rule
    """
    assert isinstance(map_dict, dict)
    if (conv_id, turn_id) not in map_dict:
        return None
    data = map_dict[(conv_id, turn_id)]
    if len(data) == 0:
        return None
    # loop over data
    matches = []
    for subrules, parsetree in data:
        if rule in subrules:
            matches.append(parsetree)
    if len(matches) == 0:
        return None
    elif len(matches) == 1:
        return matches[0]
    else:
        return matches


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
def probBoost(p_dict, r_dict, rule, method='log-odds', verbose=False):
    assert isinstance(p_dict, dict)
    assert isinstance(r_dict, dict)
    assert isinstance(rule, str)
    assert method in ['log-odds', 'diff']
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
    # print info for verbose mode
    if verbose:
        print('total pairs: {}'.format(count_pair))
        print('# in target: {}'.format(count_t))
        print('# in prime: {}'.format(count_p))
        print('# in both: {}'.format(count_pt))
    # return
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


##
# the func that counts the number of pairs that satisfy certain conditions:
# a rule's occurrance in prime is limited by prime_min and prime_max
# and its occurrance in target is limited by target_min and target_max
def count_occur(pair_dict, rule_dict, rule, prime_min, target_min, prime_max=float('inf'), target_max=float('inf')):
    """
    return: a list of tuple, in which the tuple is itself a triple tuple: [(conv_id, prime_turn_id, target_turn_id), ...]
    """
    assert isinstance(pair_dict, dict)
    assert isinstance(rule_dict, dict)
    assert isinstance(rule, str)
    assert (isinstance(prime_max, int) or isinstance(prime_max, float)) and prime_max >= 0
    assert isinstance(prime_min, int) and prime_min >= 0
    assert prime_max >= prime_min
    assert (isinstance(target_max, int) or isinstance(target_max, float)) and target_max >= 0
    assert isinstance(target_min, int) and target_min >= 0
    assert target_max >= target_min

    result = []
    for key, val in pair_dict.items():
        for pair in val:
            prime_count = rule_dict[(key, pair[0])].count(rule)
            target_count = rule_dict[(key, pair[1])].count(rule)
            if prime_count >= prime_min and prime_count <= prime_max and target_count >= target_min and target_count <= target_max:
                result.append((key, pair[0], pair[1]))
    return result


# the func that counts the number of pairs, in which a specified rule
# occurs EXACTLY ONCE in prime and target
# it is a special case of count_occur
def count_occur_once(pair_dict, rule_dict, rule, verbose=False):
    """
    find all the pairs where the specified rule occurs exactly once in both prime and target respectively
    return: [(conv_id, prime_turn_id, target_turn_id), ...]
    """
    result = count_occur(pair_dict, rule_dict, rule, prime_min=1, target_min=1, prime_max=1, target_max=1)
    # print verbose info
    if verbose:
        print(rule)
        print('# of qualified pairs: {}'.format(len(result)))
    return result


# experiment with the number of pairs, where rule appears in prime and target for once
def exp_occur():
    # load pairs_dict and rules_dict
    p_dict = pickle.load(open('swbd_dysf_adjacent_pairs.pkl', 'rb'))
    r_dict = pickle.load(open('turn_subrules_dict.pkl', 'rb'))
    # target rules
    # target_rules = ['NP -> DT NN', 'NP -> NN', 'NP -> NP PP', 'NP -> NP SBAR', 'NP -> DT', 'NP -> NNS']
    target_rules = ['NP -> DT JJ NN']
    # print info
    print('one on one repeat')
    for rule in target_rules:
        # probBoost(p_dict, r_dict, rule, verbose=True)
        res = count_occur_once(p_dict, r_dict, rule, verbose=True)
    # more than one repeat
    print('multiple on multiple repeat')
    for rule in target_rules:
        res = count_occur(p_dict, r_dict, rule, prime_min=2, target_min=2)
        print(rule)
        print('# of qualified pairs: {}'.format(len(res)))


# the func that gets the terminal nodes from a tree, filtered by the tgrep pattern
def terminal_nodes(tree_str, pattern):
    """
    tree_str: the string of a phrase structure parse tree
    pattern: the tgrep pattern
    return: a list of string of the terminal nodes (leaves)
    """
    assert isinstance(tree_str, str)
    assert isinstance(pattern, str)
    try:
        tree = ParentedTree.fromstring(tree_str)
    except Exception as e:
        print('error in constructing tree')
        raise
    else:
        res = nltk_tgrep.tgrep_nodes(tree, pattern)
        res_str = [' '.join(t.leaves()) for t in res]
        return res_str

# experiment with terminal_nodes
def exp_term_nodes():
    # print(terminal_nodes('(S (NP (DT the) (JJ big) (NN dog)) (VP bit) (NP (DT a) (NN cat)))', 'NP < (DT.NN)'))
    # NP -> NN JJ NP
    p_dict = pickle.load(open('swbd_dysf_adjacent_pairs.pkl', 'rb'))
    r_dict = pickle.load(open('turn_subrules_dict.pkl', 'rb'))
    r2t_dict = pickle.load(open('subrules_parsetree_dict.pkl', 'rb'))

    res = count_occur_once(p_dict, r_dict, 'NP -> DT JJ NN')
    parsetrees = []
    for i, item in enumerate(res):
        conv_id, prime_turn_id, target_turn_id = item
        prime_tree = ''
        target_tree = ''
        if i > 0 and prime_turn_id == res[i-1][2]:
            prime_tree = parsetrees[-1][-1]
        else:
            prime_tree = find_parsetree(r2t_dict, conv_id, prime_turn_id, 'NP -> DT JJ NN')
            target_tree = find_parsetree(r2t_dict, conv_id, target_turn_id, 'NP -> DT JJ NN')
        parsetrees.append((conv_id, prime_turn_id, target_turn_id, prime_tree, target_tree))
    #
    with open('NP<(NN.JJ.NP).csv', 'w', newline='') as fw:
        reswriter = csv.writer(fw, delimiter='\t')
        for row in parsetrees:
            reswriter.writerow(row)





# main
if __name__ == '__main__':
    # find_adj_pairs('swbd_dysf_adjacent_pairs.pkl')
    # count_rules('swbd_dysf_adjacent_pairs.pkl', 'swbd_dysf_adjacent_pairs_count.pkl')
    # reformat_count_result('swbd_dysf_adjacent_pairs_count.pkl', 'swbd_dysf_adjacent_pairs_count.txt')
    # extract_subrules_all()

    # subrules_freq('all_subrules_freq.txt')
    # turn_subrules_dict('turn_subrules_dict.pkl')
    # exp_probBoost()

    # experiment
    # subrules_parsetree_map('subrules_parsetree_dict.pkl')
    # exp_occur()
    exp_term_nodes()
