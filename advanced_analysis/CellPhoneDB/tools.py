import os 

def createDir(outdir, rm=False):
    if type(outdir).__name__ == 'list':
        for out in outdir:
            if not os.path.exists(out): 
                os.makedirs(out)
            else:
                if rm:
                    os.system(f"rm -r {out}")
                    os.makedirs(out)
    else:
        if not os.path.exists(outdir): 
            os.makedirs(outdir)
        else:
            if rm:
                os.system(f"rm -r {outdir}")
                os.makedirs(outdir)
