def main():
    n = 24
    mms = []
    bps = []
    ans = 0.0
    cnt = 0.0
    for _ in range(n):
        mm_bp = list(map(int, input().split()))
        mm = mm_bp[0]
        bp = mm_bp[1]
        mms.append(mm)
        bps.append(bp)
        ans += float(bp) / mm
        cnt += 1.0
    print(ans / cnt)
    # print(np.corrcoef(np.array(mms).flatten(), np.array(bps).flatten())[0, 1])


if __name__ == '__main__':
    main()
