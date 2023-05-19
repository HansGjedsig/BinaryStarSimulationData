def progress_bar(progress, total,start_t,t_loop):
    import time
    th=0
    tm=0
    ts=0
    if progress == 0:
        progress=1
    percent = 100*(progress/float(total))
    t1=(time.time()-start_t)+t_loop
    
    tott = t1/progress*(total-progress)
    th = int(tott/(60*60))
    tm = int((tott/(60*60)-th)*60)
    ts =  int(((tott/(60*60)-th)*60-tm)*60)
    if progress>=total:
        th, tm, ts =0,0,0
    
    print(f"\r{percent:.2f}% estimated time left: {th:02d}h {tm:02d}m {ts:02d}s",end="\r")
    return(t1)
    
def progress_bar_while(progress, total,start_t,t_loop):
    import time
    th=0
    tm=0
    ts=0
    if progress == 0:
        progress=1
    percent = 100*(progress/float(total))
    t1=(time.time()-start_t)+t_loop
    
    tott = t1/progress*(total-progress)
    th = int(tott/(60*60))
    tm = int((tott/(60*60)-th)*60)
    ts =  int(((tott/(60*60)-th)*60-tm)*60)
    if progress>=total:
        th, tm, ts =0,0,0
    
    print(f"\r{percent:.2f}% estimated time left: {th:02d}h {tm:02d}m {ts:02d}s",end="\r")
    return(t1)
