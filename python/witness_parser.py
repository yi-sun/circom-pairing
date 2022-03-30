import json

witness_file = "../build/millerloop/witness.json"
sym_file = "../build/millerloop/millerloop.sym"

sym = open(sym_file, 'r')
witness_f = open(witness_file, 'r')
witness = json.load(witness_f)

signal_list = []

for line in sym:
    count, line_num, c, signal = line.strip().split(',')
    # count is total count of signals, before optimization
    # line_num = -1 if it's optimized out I think?
    # not sure what c does
    count = int(count)
    if count == -1:
        continue
    if count < len(witness):
        signal_list.append([signal, int(witness[ count ])])

class Node(object):
    def __init__(self, value, child=None):
        self.value = value
        self.child = {}
        self.num = 0

tree_head = Node("dummy")

for signal, val in signal_list:
    signal_path = signal.split('.')
    node = tree_head
    for name in signal_path:
        if name in node.child:
            node = node.child[name]
        else:
            newnode = Node(name)
            node.child[name] = newnode
            node = newnode
    node.num = val
signal_list[0][0].split('.')

def printTree(node, level=0):
    if node != None:
        if level >= 0:
            print(' ' * 4 * level + '-> ' + node.value)
        if node.child:
            for child in node.child.values():
                printTree(child, level + 1)
        else:
            print(' ' * 4 * (level+1) + '= {}'.format(node.num))

printTree(tree_head, -1)

