# Process the results from the Java program that extract leaf nodes from the tree
# by either computing the entropy (for phrases) or simply the word frequency (for word)
# Yang Xu
# 11/23/2016

from nltk import FreqDist
import MySQLdb


# get db connection
def db_conn(db_name):
    # db init: ssh yvx5085@brain.ist.psu.edu -i ~/.ssh/id_rsa -L 1234:localhost:3306
    conn = MySQLdb.connect(host = "127.0.0.1",
                    user = "yang",
                    port = 1234,
                    passwd = "05012014",
                    db = db_name)
    return conn

# read all text of Switchboard that has dysfluencies removed
def read_swbd_text(output_file):
    conn = db_conn('swbd')
    cur = conn.cursor()

    sql = 'select rawWord from entropy_disf'
    cur.execute(sql)
    data = cur.fetchall()

    with open(output_file, 'w') as fw:
        for d in data:
            fw.write(d[0] + '\n')

# build FreqDist from a text file
def build_fdist(input_file):
    fdist = FreqDist()
    with open(input_file, 'r') as fr:
        for line in fr:
            words = line.strip().split()
            for w in words:
                fdist[w] += 1
    return fdist

##
# experiment with preliminary results
def fdist_exp():
    data_folder = '/Users/yangxu/Documents/workspace/TregexFindHeads/output/'
    data_files = ['NP -> DT NN.csv', 'NP -> NN.csv', 'NP -> NNS.csv']
    fdist = build_fdist('results/Swtichboard_dysf_rawWord.txt')

    for f in data_files:
        freqs = []
        with open(data_folder + f, 'r') as fr:
            for line in fr:
                words = line.strip().split('\t')
                freqs.append((fdist[words[0]], fdist[words[1]]))
        # write freqs to new file
        with open('results/' + f[:-4] + '_freq.txt', 'w') as fw:
            for tup in freqs:
                fw.write('\t'.join(map(str, tup)) + '\n')


####
# main
if __name__ == '__main__':
    # read_swbd_text('results/Swtichboard_dysf_rawWord.txt')
    fdist_exp()
