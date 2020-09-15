function [c, p, g, b] = gpu_checkdeadaxon(cMap1, pox, green)

  p = pox;
  g = green;
  b = 0;
  c = cMap1;
  if green <= 0
     return;
  end

  if cMap1 > 22 
      p = 0;
      g = 0;
      b = 1;
      c = 10000;
  end
end