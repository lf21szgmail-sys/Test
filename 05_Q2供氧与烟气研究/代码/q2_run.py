
import argparse
import csv
import itertools
import json
import platform
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
import PIL
from q2_load import HERE, manifest, read_book

VARS = ("Sequence", "Time", "Oxygen", "CO2", "CO")
UNITS = ("1", "s", "m3", "m3/s", "m3/s")

def quant(x):
    return np.quantile(x, [0,.1,.25,.5,.75,.9,1]).tolist() if len(x) else []

def corr(x,y):
    return float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else None

def csv_out(out, name, rows):
    with (out / (name+".csv")).open("w", encoding="utf-8-sig",newline="") as file:
        writer=csv.writer(file)
        writer.writerows(rows)

def plot(out, name, t, z, event=None):
    im=Image.new("RGB",(1000,640),"#ffffff")
    d=ImageDraw.Draw(im)
    d.text((66,10), name+" | recorded observations | t (s)", fill="#182a40")
    for j,label in enumerate(("Cumulative oxygen (m3)","CO2 (m3/s)","CO (m3/s)")):
        y=z[:,j]; top=65+j*190
        left,right=70,965
        d.rectangle((left,top,right,top+140),outline="#c1ceda")
        d.text((left,top-23),label,fill="#182a40")
        d.text((3,top),f"{max(y):.3g}",fill="#182a40")
        d.text((3,top+126),f"{min(y):.3g}",fill="#182a40")
        d.text((left,top+145),f"{t[0]:.1f}",fill="#182a40")
        d.text((right-45,top+145),f"{t[-1]:.1f}",fill="#182a40")
        if event is not None:
            xpos=left+(right-left)*(event-t[0])/(t[-1]-t[0])
            d.line((xpos,top,xpos,top+140),fill="#d36130")
        for tx,value in zip(t,y):
            x=left+(right-left)*(tx-t[0])/(t[-1]-t[0])
            p=top+138-136*(value-min(y))/(max(y)-min(y) or 1)
            d.ellipse((x-1,p-1,x+1,p+1),fill="#126998")
    im.save(out/(name+".png"))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--smoke",action="store_true")
    args=parser.parse_args()
    if args.out.resolve().is_relative_to(args.data_dir.resolve()):
        raise ValueError("Output cannot be inside source directory")
    src=list(manifest(args.data_dir))
    args.out.mkdir(parents=True,exist_ok=False)
    with (HERE/"EDA核对.csv").open(encoding="utf-8-sig",newline="") as f:
        eda={r["Heat_ID"]:r for r in csv.DictReader(f)}
    baseline=[["heat","variable","unit","n","numeric_type","blank","nonnumeric","nonfinite","first","last","min","max","zeros","zero_denominator","zero_fraction","repeats","transition_denominator","repeat_fraction","nonadjacent_duplicates_n_minus_unique","decreases","delta_quantiles_min_p10_p25_p50_p75_p90_max"]]
    intervals=[["heat","variable","kind","position","shape","excel_first","excel_last","sequence_first","sequence_last","start_s","end_s","points","point_fraction","observed_span_s","span_fraction","value","other_updates_at_entry","other_updates_at_exit","before","after"]]
    updates=[["heat","variable","observation_Hz","changing_transitions","transition_denominator","changing_fraction","observed_value_changes_per_s","inter_change_intervals_n","inter_change_interval_quantiles_s","max_equal_run_points"]]
    updates_dist=[["heat","variable","inter_change_seconds","count"]]
    sync=[["heat","oxygen_change","co2_change","co_change","count","denominator","fraction"]]
    sources=[["attachment","heat","file","sha256","rows","columns","sheet","excel_range","source_headers","eda_match"]]
    local=[["heat","event","event_excel_row","event_s","excel_row","sequence","t_s","oxygen_m3","co2_m3_s","co_m3_s","dO_m3","dCO2_m3_s","dCO_m3_s"]]
    events=[["heat","variable","rank","excel_before","excel_after","time_before_s","time_after_s","delta","before","after","other_deltas"]]
    blocks=[["heat","block","start_s","end_s","transitions","oxygen_increment","window_increment_per_s","oxygen_delta_quantiles","co2_change_count","co_change_count"]]
    relation=[["heat","gas","n","r_O_g","r_O_cumulative_g","r_dO_dg","cumulative_g_trapezoid_m3","cumulative_g_left_m3","cumulative_g_right_m3","best_positive_lag_s","r_at_lag","matched_pairs","boundary"]]
    lagrows=[["heat","gas","lag_s","pairs","r_dO_dg"]]
    for attachment,heat,path,sha,n in src:
        if args.smoke and attachment not in (1,2,5):
            continue
        labels,rows=read_book(path,heat)
        assert len(rows)==n==int(eda[heat]["Samples"])
        raw=np.array(rows,dtype=object)
        if not all(type(x) in (int,float) for x in raw.flat):
            raise ValueError(f"{heat}: blank or nonnumeric cell, audit before proceeding")
        a=raw.astype(float)
        if not np.isfinite(a).all() or not np.all(np.diff(a[:,0])==1):
            raise ValueError(f"{heat}: nonfinite cell or bad Sequence")
        t=(a[:,0]-a[0,0])*.5
        v=np.column_stack((a[:,0],t,a[:,1:]))
        d=np.diff(v,axis=0); change=d[:,2:]!=0
        duration=t[-1]
        ref=eda[heat]
        assert int(ref["Attachment"])==attachment and int(ref["Variables"])==2051
        assert np.allclose(a[[0,-1],1],[float(ref["Oxygen_First"]),float(ref["Oxygen_Last"])],rtol=0,atol=1e-9)
        for j,key in enumerate(("Oxygen","CO2","CO")):
            assert int((~change[:,j]).sum())==int(ref[key+"_Repeated_Transitions"])
            if j:
                assert int((a[:,j+1]==0).sum())==int(ref[key+"_Zero_Count"])
        sources.append([attachment,heat,path.name,sha,n,2051,heat,f"A2:BZW{n+1}",labels,True])
        print(heat,n, "EDA match",flush=True)
        for j,(key,unit) in enumerate(zip(VARS,UNITS)):
            y=v[:,j]; delta=d[:,j]
            zero=int((y==0).sum()); repeat=int((delta==0).sum())
            baseline.append([heat,key,unit,n,"integer-valued" if np.all(y==np.trunc(y)) else "float64",0,0,0,y[0],y[-1],min(y),max(y),zero,n,zero/n,repeat,n-1,repeat/(n-1),n-len(np.unique(y)),int((delta<0).sum()),quant(delta)])
            cuts=[0,*list(np.flatnonzero(delta!=0)+1),n]
            for start,stop in itertools.pairwise(cuts):
                end=stop-1
                if end==start and y[start]!=0:
                    continue
                kind="zero" if y[start]==0 else "equal"
                position="all" if start==0 and stop==n else "start" if start==0 else "end" if stop==n else "middle"
                other=[k for k in range(2,5) if k!=j]
                entry=";".join(VARS[k] for k in other if start and d[start-1,k]!=0)
                exit_=";".join(VARS[k] for k in other if stop<n and d[end,k]!=0)
                span=t[end]-t[start]; length=stop-start
                intervals.append([heat,key,kind,position,"isolated" if length==1 else "continuous",start+2,end+2,v[start,0],v[end,0],t[start],t[end],length,length/n,span,span/duration,y[start],entry,exit_,v[start-1].tolist() if start else "",v[stop].tolist() if stop<n else ""])
            gaps=np.diff(t[np.array(cuts[1:-1],dtype=int)])
            updates.append([heat,key,2,n-1-repeat,n-1,(n-1-repeat)/(n-1),(n-1-repeat)/duration,len(gaps),quant(gaps),max(np.diff(cuts))])
            for gap in np.unique(gaps):
                updates_dist.append([heat,key,gap,int((gaps==gap).sum())])
        for pattern in itertools.product((False,True),repeat=3):
            k=int(np.all(change==pattern,axis=1).sum())
            sync.append([heat,*pattern,k,n-1,k/(n-1)])
        for j in range(2,5):
            order=np.argsort(-np.abs(d[:,j]),kind="stable")[:3]
            for rank,index in enumerate(order,1):
                at=int(index+1)
                events.append([heat,VARS[j],rank,at+1,at+2,t[at-1],t[at],d[index,j],v[at-1,j],v[at,j],d[index,2:].tolist()])
                low=max(0,at-8); high=min(n,at+9)
                for row in range(low,high):
                    local.append([heat,VARS[j]+str(rank),at+2,t[at],row+2,v[row,0],t[row],*v[row,2:],*(d[row-1,2:] if row else [""]*3)])
                if rank==1 and attachment in (1,2,5,7):
                    plot(args.out,heat+"_"+VARS[j]+"_event",t[low:high],a[low:high,1:],t[at])
        if attachment==5:
            plot(args.out,heat+"_full",t,a[:,1:])
        for b,index in enumerate(np.array_split(np.arange(n-1),4),1):
            start,stop=int(index[0]),int(index[-1]+1)
            inc=a[stop,1]-a[start,1]
            blocks.append([heat,b,t[start],t[stop],len(index),inc,inc/(t[stop]-t[start]),quant(d[index,2]),int(change[index,1].sum()),int(change[index,2].sum())])
        for j in (3,4):
            y=v[:,j]
            cumul=np.r_[0,np.cumsum(.25*(y[:-1]+y[1:]))]
            lag_profiles=[]
            for lag in range(61):
                x=d[:len(d)-lag,2]; g=d[lag:,j]
                r=corr(x,g)
                lagrows.append([heat,VARS[j],lag*.5,len(x),r])
                if r is not None:
                    lag_profiles.append((lag,r,len(x)))
            lag,r,pairs=max(lag_profiles,key=lambda x:x[1])
            relation.append([heat,VARS[j],n,corr(v[:,2],y),corr(v[:,2],cumul),corr(d[:,2],d[:,j]),cumul[-1],.5*sum(y[:-1]),.5*sum(y[1:]),lag*.5,r,pairs,lag in (0,60)])
    for name,table in {"实际来源":sources,"逐炉基线":baseline,"零值与常量段":intervals,"数值更新":updates,"更新间隔分布":updates_dist,"更新同步":sync,"候选突变":events,"局部原值":local,"四段供氧":blocks,"描述关系":relation,"探索滞后曲线":lagrows}.items():
        csv_out(args.out,name,table)
    with (args.out/"运行.json").open("w",encoding="utf8") as f:
        json.dump({"python":platform.python_version(),"numpy":np.__version__,"pillow":PIL.__version__,"mode":"smoke" if args.smoke else "ten heats","baseline":"d7a5db947a3692e379f48eda6592e877803e054e","units":dict(zip(VARS,UNITS))},f,ensure_ascii=False,indent=2)

if __name__=="__main__":
    main()
