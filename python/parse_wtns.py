import argparse
import json
import pprint

def parse_signal_arrays(sigs, debug=False):
    ret = {}
    for idx, x in enumerate(sigs):
        if debug and idx % 100000 == 0:
            print(idx, len(sigs))
        var = x[0]
        val = int(x[1])
        parsed = var.split('.')
        curr = ret
        for idx, name in enumerate(parsed):
            if '[' in name:
                is_array = True
                split = name.split('[')
                array_len = len(split) - 1
                idxs = [int(x[:-1]) for x in split[1:]]
                name = split[0]
            else:
                is_array = False
                
            if name not in curr:
                if (not is_array):
                    if idx == len(parsed) - 1:
                        curr[name] = val
                    else:
                        curr[name] = {}
                        curr = curr[name]
                else:
                    if idx == len(parsed) - 1:
                        if array_len == 1:
                            curr[name] = [val]
                        elif array_len == 2:
                            curr[name] = [[val]]
                        elif array_len == 3:
                            curr[name] = [[[val]]]
                    else:
                        if sum(idxs) > 0:
                            assert()
                        else:
                            if array_len == 1:
                                curr[name] = [{}]
                                curr = curr[name][0]
                            elif array_len == 2:
                                curr[name] = [[{}]]
                                curr = curr[name][0][0]
                            elif array_len == 3:
                                curr[name] = [[[{}]]]
                                curr = curr[name][0][0][0]
            else:
                if (not is_array):
                    if idx == len(parsed) - 1:
                        assert()
                    else:
                        curr = curr[name]
                else:
                    curr_list = curr[name]
                    for dim in range(array_len - 1):
                        if len(curr_list) - 1 < idxs[dim]:
                            curr_list.append([])
                            curr_list = curr_list[-1]
                        else:
                            curr_list = curr_list[idxs[dim]]
                    if idx == len(parsed) - 1:
                        curr_list.append(val)
                    else:
                        if len(curr_list) - 1 < idxs[-1]:
                            curr_list.append({})
                            curr = curr_list[-1]
                        else:
                            curr = curr_list[idxs[-1]]
    return ret

parser = argparse.ArgumentParser()
parser.add_argument('--witness_file', type=str, default='../build/dev/witness.json')
parser.add_argument('--sym_file', type=str, default='../build/dev/dev.sym')

parser.add_argument('--debug', action='store_true', default=False)
parser.add_argument('--reparse', action='store_true', default=False)
parser.add_argument('--sig_list_file', type=str, default='../build/dev/dev.sig')

parser.add_argument('--width', type=int, default=100)
parser.add_argument('--depth', type=int, default=3)
args = parser.parse_args()

def main():
    if args.reparse:
        with open(args.witness_file, 'r') as f:
            witness = json.load(f)

        with open(args.sym_file, 'r') as f:
            signal_list = []

            for line in f:
                count, line_num, template_idx, signal = line.strip().split(',')
                # count is total count of signals, before optimization
                # line_num = -1 if it's optimized out I think?
                count = int(count)
                line_num = int(line_num)
                if line_num == -1:
                    continue
                else:
                    signal_list.append([signal, int(witness[line_num])])

        with open(args.sig_list_file, 'w') as f:
            f.write(pprint.pformat(signal_list, compact=True, width=args.width).replace("'", '"'))
    else:
        with open(args.sig_list_file, 'r') as f:
            signal_list = json.loads(f.read())

    print('Parsing')
    sigs = parse_signal_arrays(signal_list, debug=args.debug)
    pprint.pprint(sigs, compact=True, width=args.width, depth=args.depth)

if __name__ == '__main__':
    main()
