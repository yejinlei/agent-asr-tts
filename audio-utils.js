// 共享声学算法工具（纯JS，各创意页面共用）
// 包含：重采样、YIN基频、LPC+共振峰、HNR/Jitter/Shimmer、Mel谱、滤波器、窗函数
window.AU = (function(){
  function hamming(n){const w=new Float32Array(n);for(let i=0;i<n;i++)w[i]=0.54-0.46*Math.cos(2*Math.PI*i/(n-1));return w;}
  function hanning(n){const w=new Float32Array(n);for(let i=0;i<n;i++)w[i]=0.5*(1-Math.cos(2*Math.PI*i/(n-1)));return w;}
  // 重采样 linear
  function resample(x,fromSr,toSr){
    if(fromSr===toSr)return x;
    const ratio=fromSr/toSr;
    const n=Math.round(x.length/ratio);
    const y=new Float32Array(n);
    for(let i=0;i<n;i++){
      const t=i*ratio;
      const i0=Math.floor(t);const i1=Math.min(i0+1,x.length-1);
      const frac=t-i0;
      y[i]=x[i0]*(1-frac)+x[i1]*frac;
    }
    return y;
  }
  // YIN 基频估计
  function yin(frame,sr,threshold=0.1,fmin=60,fmax=500){
    const N=frame.length;
    const buf=new Float32Array(N/2);
    // diff
    for(let tau=0;tau<N/2;tau++){
      let s=0;
      for(let i=0;i<N/2;i++){const d=frame[i]-frame[i+tau];s+=d*d;}
      buf[tau]=s;
    }
    // cmndf
    buf[0]=1;let acc=0;
    for(let tau=1;tau<N/2;tau++){acc+=buf[tau];buf[tau]*=tau/acc;}
    // find minimum below threshold
    const tauMin=Math.floor(sr/fmax);const tauMax=Math.floor(sr/fmin);
    let minTau=-1,minVal=1;
    for(let tau=tauMin;tau<Math.min(tauMax,N/2);tau++){
      if(buf[tau]<threshold&&buf[tau]<minVal){minTau=tau;minVal=buf[tau];}
    }
    if(minTau<0)return {f0:0,voiced:false};
    // parabola
    let f=minTau;
    if(f>0&&f<N/2-1){
      const s0=buf[f-1],s1=buf[f],s2=buf[f+1];
      const d=(s2-s0)/(2*(2*s1-s2-s0));
      f=f+d;
    }
    return {f0:sr/f,voiced:true};
  }
  // 自相关
  function autocorr(x){const N=x.length;const r=new Float32Array(N);
    for(let k=0;k<N;k++){let s=0;for(let i=0;i<N-k;i++)s+=x[i]*x[i+k];r[k]=s;}
    return r;
  }
  // Levinson-Durbin
  function levinson(r,p){
    const a=new Float32Array(p+1);const k=new Float32Array(p);
    a[0]=1;let E=r[0];
    for(let i=1;i<=p;i++){
      let lambda=-r[i];
      for(let j=1;j<i;j++)lambda-=a[j]*r[i-j];
      lambda/=E;
      k[i-1]=-lambda;a[i]=lambda;
      for(let j=1;j<Math.floor(i/2)+1;j++){
        const aj=a[j],ajp=a[i-j];
        a[j]=aj+lambda*ajp;
        a[i-j]=ajp+lambda*aj;
      }
      E*=(1-lambda*lambda);
    }
    return{a,k,E};
  }
  // 解多项式根（伴随矩阵特征值近似，使用Durand-Kerner）
  function polyRoots(a){
    const n=a.length-1;
    // 单根近似：找共振峰只用极角
    const roots=[];
    // Durand-Kerner
    const z=[];for(let i=0;i<n;i++)z.push({r:Math.cos(2*Math.PI*i/n),i:Math.sin(2*Math.PI*i/n)});
    // 归一化
    const a0=a[0];for(let i=0;i<a.length;i++)a[i]/=a0;
    for(let iter=0;iter<80;iter++){
      let md=0;
      for(let i=0;i<n;i++){
        let p={r:1,i:0};
        for(let k=0;k<n;k++){// p = a[0]*z^n+...+a[n]
          // multiply by z[i]
          const nr=p.r*z[i].r-p.i*z[i].i,ni=p.r*z[i].i+p.i*z[i].r;
          p.r=nr;p.i=ni;p.r+=a[k+1];
        }
        // product (z[i]-z[j])
        let d={r:1,i:0};
        for(let j=0;j<n;j++)if(j!==i){
          const nr=d.r*(z[i].r-z[j].r)-d.i*(z[i].i-z[j].i);
          const ni=d.r*(z[i].i-z[j].i)+d.i*(z[i].r-z[j].r);
          d.r=nr;d.i=ni;
        }
        // divide
        const den=d.r*d.r+d.i*d.i;
        const dv=(p.r*d.r+p.i*d.i)/den,di=(p.i*d.r-p.r*d.i)/den;
        z[i].r-=dv;z[i].i-=di;md=Math.max(md,Math.sqrt(dv*dv+di*di));
      }
      if(md<1e-6)break;
    }
    for(let i=0;i<n;i++){
      const re=z[i].r,im=z[i].i;
      const mag=Math.sqrt(re*re+im*im);
      if(mag>0.6&&mag<1.05&&im>0.01){
        const ang=Math.atan2(im,re);
        roots.push({mag,ang});
      }
    }
    return roots;
  }
  // LPC 共振峰：帧→F1,F2,F3,F4
  function formants(frame,sr,p=12){
    // preemph
    const N=frame.length;
    const pre=new Float32Array(N);
    pre[0]=frame[0];for(let i=1;i<N;i++)pre[i]=frame[i]-0.67*frame[i-1];
    // window
    const w=hamming(N);
    const wf=new Float32Array(N);for(let i=0;i<N;i++)wf[i]=pre[i]*w[i];
    // autocorr
    const r=autocorr(wf);
    // levinson
    const{ a }=levinson(r.slice(0,p+1),p);
    // roots
    const roots=polyRoots(a);
    const F=roots.map(rt=>{
      const f=rt.ang*sr/(2*Math.PI);
      const bw=-0.5*(sr/(2*Math.PI))*Math.log(rt.mag);
      return{f,bw};
    }).filter(x=>x.f>50&&x.f<sr/2-100&&x.bw<500);
    F.sort((a,b)=>a.f-b.f);
    return F;
  }
  // RMS energy
  function rms(x){let s=0;for(let i=0;i<x.length;i++)s+=x[i]*x[i];return Math.sqrt(s/x.length);}
  // HNR（自相关法简化）
  function hnr(frame,sr){
    const N=frame.length;
    const w=hamming(N);const wf=new Float32Array(N);for(let i=0;i<N;i++)wf[i]=frame[i]*w[i];
    const r=autocorr(wf);
    // R0 = r[0] energy; find max in [2ms..20ms]
    const tmin=Math.floor(sr*0.002),tmax=Math.floor(sr*0.02);
    let mx=0,mxi=tmin;
    for(let i=tmin;i<Math.min(tmax,N);i++)if(r[i]>mx){mx=r[i];mxi=i;}
    if(r[0]<1e-12)return 0;
    const ratio=mx/(r[0]-mx);
    if(ratio<=0)return -20;
    return 10*Math.log10(ratio);
  }
  // 从AudioBuffer得到单声道float32
  function toMono(buf){
    let x;
    if(buf.numberOfChannels===1){x=buf.getChannelData(0).slice();}
    else{
      const L=buf.getChannelData(0),R=buf.getChannelData(1);
      x=new Float32Array(L.length);
      for(let i=0;i<L.length;i++)x[i]=(L[i]+R[i])/2;
    }
    return x;
  }
  // 分析文件或buffer：返回 metrics
  async function analyzeBuffer(buf,sr=16000){
    const x=toMono(buf);
    const x16=resample(x,buf.sampleRate,sr);
    const frameSize=Math.floor(sr*0.04),hop=Math.floor(sr*0.01);
    const frames=[];
    let totalRMS=0,voicedCount=0;
    const F0s=[],F1s=[],F2s=[],HNRs=[];
    for(let i=0;i+frameSize<x16.length;i+=hop){
      const f=x16.subarray(i,i+frameSize);
      const e=rms(f);totalRMS+=e;
      const y=yin(f,sr,0.12,70,500);
      const F=formants(f,sr,12);
      const h=hnr(f,sr);
      frames.push({rms:e,f0:y.f0,voiced:y.voiced,F,HNR:h});
      if(y.voiced){voicedCount++;F0s.push(y.f0);HNRs.push(h);}
      if(F[0])F1s.push(F[0].f);
      if(F[1])F2s.push(F[1].f);
    }
    const avg=arr=>arr.length?arr.reduce((a,b)=>a+b,0)/arr.length:0;
    const med=arr=>{if(!arr.length)return 0;const s=arr.slice().sort((a,b)=>a-b);return s[Math.floor(s.length/2)];};
    return{
      sr,totalSamples:x16.length,duration:x16.length/sr,
      avgRMS:totalRMS/Math.max(1,frames.length),
      voicedRatio:voicedCount/Math.max(1,frames.length),
      avgF0:med(F0s),avgF1:avg(F1s),avgF2:avg(F2s),avgHNR:avg(HNRs),
      frames,raw:x16
    };
  }
  // Jitter & Shimmer（对浊音帧）
  function jitterShimmer(frames){
    const voiced=frames.filter(f=>f.voiced&&f.f0>50&&f.f0<500);
    if(voiced.length<3)return{jitter:0,shimmer:0};
    let js=0,ss=0;const n=voiced.length;
    for(let i=1;i<n;i++){
      const T0=1/voiced[i-1].f0,T1=1/voiced[i].f0;
      js+=Math.abs(T1-T0)/((T0+T1)/2);
      ss+=Math.abs(voiced[i].rms-voiced[i-1].rms)/((voiced[i-1].rms+voiced[i].rms)/2+1e-9);
    }
    return{jitter:js/(n-1)*100,shimmer:ss/(n-1)*100};
  }
  // 生成正弦波
  function sine(freq,sr,dur,phase=0){
    const n=Math.floor(sr*dur);
    const x=new Float32Array(n);
    for(let i=0;i<n;i++)x[i]=Math.sin(2*Math.PI*freq*i/sr+phase);
    return x;
  }
  // 音频播放
  function playFloat(x,sr){
    const ctx=new(window.AudioContext||window.webkitAudioContext)();
    const b=ctx.createBuffer(1,x.length,sr);
    b.getChannelData(0).set(x);
    const src=ctx.createBufferSource();src.buffer=b;src.connect(ctx.destination);src.start();
    return src;
  }
  // 简单FIR低通(窗函数法)
  function designLPF(cutoff,sr,N=101){
    const fc=cutoff/sr;
    const h=new Float32Array(N);
    const w=hamming(N);
    const mid=Math.floor(N/2);
    for(let n=0;n<N;n++){
      const k=n-mid;
      h[n]=k===0?2*fc:Math.sin(2*Math.PI*fc*k)/(Math.PI*k);
      h[n]*=w[n];
    }
    return h;
  }
  function applyFIR(x,h){
    const y=new Float32Array(x.length);
    for(let i=0;i<x.length;i++){
      let s=0;
      for(let k=0;k<h.length;k++){if(i-k>=0)s+=h[k]*x[i-k];}
      y[i]=s;
    }
    return y;
  }
  // FFT (Cooley-Tukey radix-2)
  function fft(re,im){
    const n=re.length;
    // bitreverse
    for(let i=1,j=0;i<n;i++){
      let bit=n>>1;
      for(;j&bit;bit>>=1)j^=bit;
      j^=bit;
      if(i<j){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}
    }
    for(let len=2;len<=n;len<<=1){
      const ang=-2*Math.PI/len;
      const wlr=Math.cos(ang),wli=Math.sin(ang);
      for(let i=0;i<n;i+=len){
        let wr=1,wi=0;
        for(let k=0;k<len/2;k++){
          const ur=re[i+k],ui=im[i+k];
          const vr=re[i+k+len/2]*wr-im[i+k+len/2]*wi;
          const vi=re[i+k+len/2]*wi+im[i+k+len/2]*wr;
          re[i+k]=ur+vr;im[i+k]=ui+vi;
          re[i+k+len/2]=ur-vr;im[i+k+len/2]=ui-vi;
          const nwr=wr*wlr-wi*wli;nwi=wr*wli+wi*wlr;wr=nwr;wi=nwi;
        }
      }
    }
  }
  // 下一个2的幂
  function nextPow2(n){let p=1;while(p<n)p<<=1;return p;}
  // 绘制波形到canvas
  function drawWave(canvas,x,color='#60a5fa',bg='#0f172a'){
    const ctx=canvas.getContext('2d');
    const W=canvas.width,H=canvas.height;
    ctx.fillStyle=bg;ctx.fillRect(0,0,W,H);
    ctx.strokeStyle=color;ctx.lineWidth=1;
    ctx.beginPath();
    const step=Math.ceil(x.length/W);
    for(let i=0;i<W;i++){
      const st=i*step,en=Math.min(st+step,x.length);
      let mx=-1e9,mn=1e9;
      for(let j=st;j<en;j++){if(x[j]>mx)mx=x[j];if(x[j]<mn)mn=x[j];}
      const y1=H/2-(mx*H*0.45);
      const y2=H/2-(mn*H*0.45);
      if(i===0)ctx.moveTo(i,y1);else ctx.lineTo(i,y1);
      ctx.lineTo(i,y2);
    }
    ctx.stroke();
    // center line
    ctx.strokeStyle='rgba(255,255,255,0.15)';ctx.beginPath();ctx.moveTo(0,H/2);ctx.lineTo(W,H/2);ctx.stroke();
  }
  // 绘制频谱图
  function drawSpectrum(canvas,re,im,N,sr,maxF=8000){
    const ctx=canvas.getContext('2d');
    const W=canvas.width,H=canvas.height;
    ctx.fillStyle='#0f172a';ctx.fillRect(0,0,W,H);
    const bins=Math.min(N/2,Math.floor(maxF/sr*N));
    // find max
    let mx=-1e9;
    for(let i=0;i<bins;i++){const m=Math.sqrt(re[i]*re[i]+im[i]*im[i]);if(m>mx)mx=m;}
    ctx.strokeStyle='#60a5fa';ctx.lineWidth=1.5;ctx.beginPath();
    for(let i=0;i<W;i++){
      const b=Math.floor(i/W*bins);
      const m=Math.sqrt(re[b]*re[b]+im[b]*im[b]);
      const db=20*Math.log10(m/(mx+1e-9)+1e-9);
      const y=H/2-(db+60)/60*H/2;
      if(i===0)ctx.moveTo(i,y);else ctx.lineTo(i,y);
    }
    ctx.stroke();
  }
  // 文件到AudioBuffer
  async function decodeFile(file){
    const arr=await file.arrayBuffer();
    const ctx=new(window.AudioContext||window.webkitAudioContext)();
    return await ctx.decodeAudioData(arr);
  }
  // mic
  async function startMic(onData){
    const stream=await navigator.mediaDevices.getUserMedia({audio:{echoCancellation:false,noiseSuppression:false}});
    const ctx=new(window.AudioContext||window.webkitAudioContext)({sampleRate:16000});
    const src=ctx.createMediaStreamSource(stream);
    const proc=ctx.createScriptProcessor(1024,1,1);
    proc.onaudioprocess=e=>onData(e.inputBuffer.getChannelData(0));
    src.connect(proc);proc.connect(ctx.destination);
    return{stop(){stream.getTracks().forEach(t=>t.stop());proc.disconnect();ctx.close();},ctx};
  }
  return{hamming,hanning,resample,yin,autocorr,levinson,polyRoots,formants,rms,hnr,toMono,analyzeBuffer,jitterShimmer,sine,playFloat,designLPF,applyFIR,fft,nextPow2,drawWave,drawSpectrum,decodeFile,startMic};
})();
