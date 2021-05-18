function cmap = redblue(n)

vlin = linspace(0,1,n/2)'; vinv = flipud(vlin); vone = 1+0*vlin;
cmap = [vlin,vlin,vone ;
        vone,vinv,vinv ];