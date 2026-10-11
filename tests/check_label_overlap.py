"""Detecta solapes reales: para cada figura, vuelve a dibujarla y compara las cajas
de todos los objetos de texto dentro de cada eje."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, sys
def report(fig, name):
    fig.canvas.draw(); r=fig.canvas.get_renderer(); bad=[]
    items=[]
    for ax in fig.axes:
        vis_axis = ax.axison and ax.get_visible()
        cand = list(ax.texts) + [ax.title]
        if vis_axis:
            cand += list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + [ax.xaxis.label, ax.yaxis.label]
        for t in cand:
            if not t.get_text().strip() or not t.get_visible(): continue
            try: bb=t.get_window_extent(r)
            except Exception: continue
            items.append((ax,t,bb))
        leg=ax.get_legend()
        if leg is not None:
            try: items.append((ax,leg,leg.get_window_extent(r)))
            except Exception: pass
    for i in range(len(items)):
        for j in range(i+1,len(items)):
            a,ta,ba=items[i]; b,tb,bb=items[j]
            if a is not b: continue
            ov=(min(ba.x1,bb.x1)-max(ba.x0,bb.x0))*(min(ba.y1,bb.y1)-max(ba.y0,bb.y0))
            if ov>8 and min(ba.x1,bb.x1)>max(ba.x0,bb.x0) and min(ba.y1,bb.y1)>max(ba.y0,bb.y0):
                sa=ta.get_text() if hasattr(ta,"get_text") else "legend"
                sb=tb.get_text() if hasattr(tb,"get_text") else "legend"
                bad.append((sa[:28],sb[:28],int(ov)))
    if bad:
        print(f"\n--- {name}: {len(bad)} solapes")
        for s in bad[:12]: print(f"    '{s[0]}'  x  '{s[1]}'   area {s[2]}")
    return len(bad)
