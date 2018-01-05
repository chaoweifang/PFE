function seg = Ucm2Seg(ucm2,thre)
labels2 = bwlabel(ucm2 <= thre);
seg = labels2(2:2:end, 2:2:end);
end