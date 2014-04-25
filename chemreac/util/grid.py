def bounds(nstencil, N):
    # appropriate when lrefl=False, refl=False
    return [(
        max(0, min(N-nstencil, i - (nstencil-1)//2)),
        min(N, max(  nstencil, i + (nstencil+1)//2))
    ) for i in range(N)]
