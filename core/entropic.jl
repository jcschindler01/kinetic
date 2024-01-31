



module KEntropic

using DataStructures

## histogram with even bin widths
function hist(x; nbins=10, xmin=0, xmax=1)
	h = zeros(Int, nbins)
	dx = (xmax-xmin)/nbins
	for i=1:length(x)
		h[1+floor(Int,x[i]/dx)]+=1
	end
	return h
end

## 2d histogram
function hist2d(x, y; nbins=3, xymin=0, xymax=1)
	h = zeros(Int, (nbins,nbins))
	dxy = (xymax-xymin)/nbins
	for i=1:length(x)
		h[1+floor(Int,x[i]/dxy), 1+floor(Int,y[i]/dxy)]+=1
	end
	return h
end

## return qref
function qref(qdict, macrostring)
	if macrostring in keys(qdict)
		return qdict[macrostring]/qdict["_N"]
	else
		return 1/qdict["_N"]
	end
end





end




