// Shared JS: TOC highlight + Tab switching (auto-init on load)
(function(){
  function init(){
    const toc=document.querySelector('aside.toc nav');
    if(toc){
      const links=[...toc.querySelectorAll('a')];
      const secs=links.map(a=>{const s=document.querySelector(a.getAttribute('href'));return s;}).filter(Boolean);
      const primary=getComputedStyle(document.documentElement).getPropertyValue('--primary').trim()||'#0ea5e9';
      const light=getComputedStyle(document.documentElement).getPropertyValue('--primary-l').trim()||'#e0f2fe';
      function upd(){
        const y=window.scrollY+120;let cur=secs[0];
        for(const s of secs){if(s&&s.offsetTop<=y)cur=s;}
        links.forEach(a=>{
          const on=a.getAttribute('href')==='#'+(cur&&cur.id);
          a.classList.toggle('on',on);
          if(on){a.style.color=primary;a.style.borderLeftColor=primary;a.style.background=light;a.style.fontWeight='600';}
          else{a.style.color='';a.style.borderLeftColor='';a.style.background='';a.style.fontWeight='';}
        });
      }
      window.addEventListener('scroll',upd,{passive:true});upd();
    }
    // tabs
    document.querySelectorAll('.tabs').forEach(tg=>{
      const panels=[];
      let sib=tg.nextElementSibling;
      while(sib && sib.classList.contains('tab-pnl')){panels.push(sib);sib=sib.nextElementSibling;}
      function activate(k){
        tg.querySelectorAll('.tab').forEach(b=>{
          const on=b.dataset.t===k;
          b.classList.toggle('on',on);
          b.setAttribute('aria-selected',on?'true':'false');
        });
        panels.forEach(pn=>pn.classList.toggle('on',pn.dataset.p===k));
      }
      tg.querySelectorAll('.tab').forEach(btn=>{
        btn.addEventListener('click',()=>activate(btn.dataset.t));
      });
      // set initial state based on first .on tab
      const initial=tg.querySelector('.tab.on');
      if(initial)activate(initial.dataset.t);
    });
  }
  if(document.readyState==='loading'){document.addEventListener('DOMContentLoaded',init);}
  else{init();}
})();
