import sys

EPSILON = 1e-5

def main():
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <file1> <file2> [debug]')
        return 1
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    debug = len(sys.argv) > 3
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        joint_iter = zip(f1, f2)
        next(joint_iter)  # skip execution time
        for i, (l1, l2) in enumerate(joint_iter):
            if l1 and l2:
                x1, y1 = l1.split(',')
                x2, y2 = l2.split(',')
                if abs(float(x1) - float(x2)) > EPSILON or abs(float(y1) - float(y2)) > EPSILON:
                    if debug:
                        print(f'Planet {i} differs')
                    return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
