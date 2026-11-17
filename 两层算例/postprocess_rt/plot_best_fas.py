"""绘制已保存最优组合的加速度FAS；不重新选站或覆盖反透射系数结果。"""
from pathlib import Path
import argparse
import copy
import csv
import json

import numpy as np
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import postprocess_rt as rt


def fas(x, dt, nfft, taper_fraction):
    # Positive-frequency Fourier-transform magnitude, not PSD or 2|FFT|/N.
    return dt*np.abs(np.fft.rfft(x*tukey(len(x), alpha=taper_fraction), n=nfft))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=Path(__file__).resolve().parents[1]/'算例'/'Planewave_upper_15_P')
    parser.add_argument('--max-frequency',type=float,default=300.)
    args=parser.parse_args();case=args.case.resolve();process=case/'result'/'过程文件'
    if not np.isfinite(args.max_frequency) or args.max_frequency<=0:raise ValueError('max-frequency must be positive')
    report=json.loads((process/'run_report.json').read_text(encoding='utf-8-sig'))
    saved=json.loads((process/'processing_config.json').read_text(encoding='utf-8-sig'))
    o=rt.settings();o.update(saved)
    if o['trace_quantity']!='acceleration':raise ValueError('本脚本输出加速度FAS，需要acceleration记录。')
    model=rt.load_model(case/'prepare');index=rt.TraceIndex(case/'sema',o)
    best=[p for p in report['best_station_pairs'] if p['best_status']=='selected']
    if not best:raise ValueError('保存结果中没有最优台站组合。')
    lookup={(r['family'],r['phase'],r['station']):r for r in rt.read_csv(process/'amplitudes.csv')}
    raw=copy.deepcopy(o);raw.update(bandpass_hz=None,detrend='none',baseline_interval_s=None)
    records=[];dt=None;nfft=None;frequency=None
    for pair in best:
        for role,phase,name in [('incident','I_'+model.config['incident_wave'],pair['incident_station']),('outgoing',pair['phase'],pair['outgoing_station'])]:
            row=lookup[pair['family'],phase,name]
            t,signal,_,_=rt.project_trace(index,name,model.phases[phase],raw)
            tf,filtered,_,_=rt.project_trace(index,name,model.phases[phase],o)
            np.testing.assert_array_equal(t,tf)
            current_dt=float(t[1]-t[0])
            if dt is None:
                dt=current_dt;nfft=2**int(np.ceil(np.log2(len(t)*8)))
                frequency=np.fft.rfftfreq(nfft,dt)
            if not np.isclose(dt,current_dt,rtol=1e-9):raise ValueError('不同台站采样间隔不一致。')
            if len(t)>nfft:raise ValueError('完整记录长度超过公共FFT长度。')
            start,end=float(row['start_s']),float(row['end_s']);mask=(t>=start)&(t<=end)
            metrics,_,_=rt.measure(t,filtered,start,end,o)
            for key in ('direct_peak','hilbert_peak'):
                np.testing.assert_allclose(metrics[key],float(row[key]),rtol=1e-9)
                np.testing.assert_allclose(metrics[key],float(pair[role+'_'+key]),rtol=1e-9)
            records.append(dict(pair_phase=pair['phase'],role=role,target_phase=phase,station=name,
                                start_s=start,end_s=end,full_start_s=float(t[0]),full_end_s=float(t[-1]),
                                samples=len(t),window_samples=int(mask.sum()),dt=dt,
                                full_raw=fas(signal,dt,nfft,o['taper_fraction']),
                                window_raw=fas(signal[mask],dt,nfft,o['taper_fraction']),
                                window_filtered=fas(filtered[mask],dt,nfft,o['taper_fraction'])))
    out=case/'result'/'最优台站FAS';out.mkdir(parents=True,exist_ok=True)
    rt.configure_plots()
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    colors={'incident':'#246a9a','outgoing':'#cf6531'}
    labels={'R_P':'反射 P','T_P':'透射 P','T_SV':'透射 SV','R_SV':'反射 SV'}
    maximum=min(args.max_frequency,frequency[-1]);display=frequency<=maximum
    def panel(ax,pair,key):
        group=[r for r in records if r['pair_phase']==pair['phase']]
        for r in group:
            label=('入射' if r['role']=='incident' else '出射')+'：'+r['station']
            if key!='full_raw':label+=f"\n窗 [{r['start_s']:.4f}, {r['end_s']:.4f}] s"
            ax.plot(frequency[display],r[key][display],lw=1.7,color=colors[r['role']],label=label)
        ax.set(xlim=(0,maximum),ylim=(0,None),xlabel='频率 / Hz',ylabel='加速度 FAS / (m/s)')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
        ax.axvline(model.config['frequency_hz'],color='#777',ls=':',lw=1)
        if o['bandpass_hz']:
            ax.axvspan(*o['bandpass_hz'],color='#c7d6df',alpha=.25,label='原处理带通频带')
        ax.legend(loc='upper right',fontsize=8.5);ax.grid(alpha=.16)
    def save(fig,name):
        for ext in ('png','pdf'):fig.savefig(out/(name+'.'+ext),dpi=190,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    fig,axes=plt.subplots(len(best),1,figsize=(14,3.6*len(best)+1),layout='constrained',squeeze=False)
    for ax,pair in zip(axes[:,0],best):
        panel(ax,pair,'full_raw')
        ax.set_title(f"{labels.get(pair['phase'],pair['phase'])}组合 | 保存的Hilbert系数={float(pair['hilbert_ratio_abs']):.6f}，误差={pair['best_error_pct']:.4f}%",loc='left',fontsize=12)
    fig.suptitle(f"{case.name}｜三组最优台站：全记录未滤波 FAS\n按各自目标P/SV极化方向投影；两端Tukey渐消，不逐道归一化",fontsize=16)
    fig.supxlabel('全记录包含所有到达波；阴影仅标示原90–110 Hz处理频带，不代表本图已滤波。' if o['bandpass_hz']==[90.,110.] else '全记录包含所有到达波；阴影仅标示原处理频带，不代表本图已滤波。',fontsize=10)
    save(fig,'最优三组台站_FAS_全记录')
    fig,axes=plt.subplots(len(best),2,figsize=(17,3.8*len(best)+1),layout='constrained',squeeze=False)
    for axes_row,pair in zip(axes,best):
        for ax,key,title in zip(axes_row,['window_raw','window_filtered'],['目标拾取窗：未滤波','目标拾取窗：原处理参数']):
            panel(ax,pair,key);ax.set_title(labels.get(pair['phase'],pair['phase'])+'｜'+title,fontsize=12,loc='left')
    fig.suptitle(f"{case.name}｜最优组合拾取窗 FAS 对照\n沿用已保存的实际拾取窗；右列先对全记录按原设置预处理，再截窗和渐消",fontsize=15)
    fig.supxlabel('FAS = Δt × |FFT(w·a)|；只展示正频率，不乘2，不作相干增益补偿；补零只加密曲线，不提高实际频率分辨率。',fontsize=10)
    save(fig,'最优三组台站_FAS_拾取窗对照')
    for r in records:
        with (out/(r['pair_phase']+'_'+r['role']+'_FAS.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.writer(f);writer.writerow(['frequency_Hz','full_raw_FAS_m_per_s','window_raw_FAS_m_per_s','window_processed_FAS_m_per_s'])
            writer.writerows(zip(frequency,r['full_raw'],r['window_raw'],r['window_filtered']))
    metadata=dict(case=str(case),selection_source=str(process/'run_report.json'),best_pairs=best,
                  bandpass_hz=o['bandpass_hz'],filter_order=o['filter_order'],taper_fraction=o['taper_fraction'],
                  definition='dt*abs(rfft(tukey(N,alpha)*acceleration)); positive frequencies, no doubling, no normalization or taper gain correction',
                  nfft=nfft,frequency_bin_spacing_hz=float(frequency[1]),plot_max_frequency_hz=maximum,
                  records=[{k:v for k,v in r.items() if not isinstance(v,np.ndarray)} for r in records],
                  verification='Reproduced all six direct and Hilbert peaks against saved amplitude rows and best-pair records.')
    (out/'FAS计算说明.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Saved:',out)
    print('Verified six saved amplitude picks; full-spectrum CSV range:',frequency[0],frequency[-1])


if __name__=='__main__':main()
