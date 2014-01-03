   
function [depth, medianCurveId] = wfbplot(data, x, weight, params)
% Produces functional boxplots or enhanced functional boxplots of the given functional data. 
% It can also be used to carry out functional data ordering based on band depth. 
%
%Params
%
%	data: a p-by-n functional data matrix where n is the number of curves, and p is defined below.
%         alternatively, a functional data object or functional parameter
%         object can also be specified. 
%	x:	the x coordinates of curves. Defaults to 1:p where p is the number of x coordinates.
%	method: the method to be used to compute band depth. Can be one of "BD2", "BD3", "MBD" or "Both" 
%	        with a default of "MBD".
%	depth:	a vector giving band depths of curves. If missing, band depth computation is conducted.
%	showOut: logical. If TRUE (the default) then the outliers will be shown. 
%	centerColor: color of the central region. Defaults to be magenta.
%	barColor: color of bars in a functional boxplot. Defaults to be blue.
%	factor: the constant factor to inflate the middle box and determine fences for outliers. Defaults to 
%	        be 1.5 as in a classical boxplot.
%
%Details
%
%	For functional data, the band depth (BD) or modifed band depth (MBD) allows for ordering a sample of 
%	curves from the center outwards and, thus, introduces a measure to define functional quantiles and 
%	the centrality or outlyingness of an observation. A smaller rank is associated with a more central 
%	position with respect to the sample curves. "BD2" uses two curves to determine a band and "BD3" uses 
%	three curves. BD usually provides many ties (curves have the same depth values), but MBD does not. 
%	The method "Both" uses BD2 first and then uses MBD to break ties.
%
%Value
%
%	depth:	band depths of given curves.
%	medianCurveId: column indices of detected outliers.
%
%Functional boxplot's Author(s)
%
%	Ying Sun: sunwards@stat.tamu.edu
%
%	Marc G. Genton: genton@stat.tamu.edu
%
%Weighted functional boxplot's Author(s)
%
%   Yi Hong, Brad Davis, J.S. Marron, Roland Kwitt, Marc Niethammer
%   Yi Hong: yihong@cs.unc.edu
%
%References
%
%	Sun, Y. and Genton, M. G. (2011), "Functional Boxplots," Journal of Computational and Graphical Statistics,
%	to appear.
%
%	Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for functional data," Journal of the American
%	Statistical Association, 104, 718-734.
%
%   Y. Hong, B. Davis, J.S. Marron, R. Kwitt, and M. Niethammer (2013),
%   "Weighted functional boxplots with application to statistical atlas
%   construction", MICCAI 2013.
%

%%default values
if (~isfield(params, 'factor'))       factor = 1.5;       end
if (~isfield(params, 'barColor'))     barcol = 'b';       end
if (~isfield(params, 'centerColor'))  color = 'm';        end       
if (~isfield(params, 'showOut'))      showOut = 'True';   end
if (~isfield(params, 'method'))       method = 'MBD';     end
if (~isfield(params, 'depth'))        depth = [];         end

% If data is an fd object or fdPar object extract it
if(isa_fdPar(data)), data = getfd(data); end
if(isa_fd(data))
   if isempty(x) || length(x)==1,
      rr = getbasisrange(getbasis(data));
      if isempty(x), npt = 101;
      else npt = x; end
      x = linspace(rr(1),rr(2),npt)';
   end
   data = eval_fd(x,data);
end

% default value for x
[tp,n]=size(data);
if isempty(x), x = (1:tp)'; end
if (length(x) ~= tp), error('Dimensions of data and x do not match'); end

%compute band depth	
if isempty(depth)
   if strcmp(params.method,'BD2')
      depth=BD2(data', weight);
   elseif strcmp(params.method,'BD3')
      depth=BD3(data', weight);
   elseif strcmp(params.method,'MBD')
      depth=MBD(data', weight);
   elseif strcmp(params.method,'Both')
      depth=round(BD2(data', weight)*10000)+MBD(data', weight);  
   end
end
    
[dp_s,index]=sort(depth,'descend');

% compute the 50% center region and the 99.3% fences based on weights
weight_sum = 0;
index_50_flag = 0;
index_993_flag = 0;
for iI = 1:length(index)
    weight_sum = weight_sum + weight( index( iI ) );
    if weight_sum >= 0.993 && index_993_flag == 0
        index_993 = iI;
        index_993_flag = 1;
	elseif weight_sum >= 0.5 && index_50_flag == 0
        index_50 = iI;
        index_50_flag = 1;
    end	
end
m = index_50;
center = data(:,index(1:m));
inf = min(center,[],2)';
sup = max(center,[],2)';

% compute outliers
dist = params.factor*(sup-inf);
upper = sup+dist;
lower = inf-dist;      
outly=sum(or(data<lower'*ones(1,n),data>upper'*ones(1,n)));  % 1.5IQR
outly(index(index_993+1:end)) = 1;  % 99.3%
outpoint = find(outly);
out = data(:,outpoint);
weight_out = weight(outpoint);
good = data;
good(:,outpoint) = [];
maxcurve=max(good,[],2)';
mincurve=min(good,[],2)';

% plot the outliers
if params.showOut
    if sum(outly) > 0
        hold all;
        weight_out = weight_out ./ sum( weight_out);
        for iOut = 1:length(weight_out)
            plot( x, out(:, iOut), 'Color', [1.0-weight_out(iOut) 1.0-weight_out(iOut) ...
                  1.0-weight_out(iOut)], 'LineStyle', '--' );
        end
    end
end

% compute the position of the bars
barval = (x(1)+x(tp))/2;
loc = find(sort([x;barval])==barval);
bar = loc(1);

% plot the 
hold all;
[xinv,xindex] = sort(x,'descend');
xx = [x;xinv];
supinv = sup(xindex);
yy = [inf,supinv];
h = fill(xx,yy,params.centerColor);
set(h,'edgecolor',params.barColor,'LineWidth',2);
hold all;
line([x(bar) x(bar)],[maxcurve(bar) sup(bar)],'Color',params.barColor,'LineWidth',2);
hold all;
line([x(bar) x(bar)],[mincurve(bar) inf(bar)],'Color',params.barColor,'LineWidth',2);
       
% plot the median, minimal, maximum, curve    
hold all;
plot(x,data(:,index(1)),'Color','k','LineWidth',2);
plot(x,maxcurve,'Color','b','LineWidth',2);
plot(x,mincurve,'Color','b','LineWidth',2);
hold off;

medianCurveId = index(1);
end


function combinat = combinat(n,p)
if n<p 
    combinat=0;
else
    combinat=nchoosek(n,p);
end
end

function resultados = estaEntre(v,matrizDatos)
% The input in this function is the data matrix and the pair of indexes and
% the output is a vector of dimension n, where the coordinate k takes the 
% value 1 or 0 whether the observation k is inside the band or not.
[n,p]=size(matrizDatos);                      
Z=matrizDatos;
inf=min(Z(v,:))';
sup=max(Z(v,:))';
resultados=sum((and(Z'<=sup*ones(1,n),Z'>=inf*ones(1,n))))==p;
end

function [contg]=BD2(matrizDatos, weight)
% With this function we calculate the band depth of every observation in 
% the matrix of data: matrizDatos. The rows are the observations and the 
% columns are the variables.
[n,p]=size(matrizDatos);
cont=zeros(1,n); % The depths are initialized to 0.
sum_weight = 0;
for i=1:(n-1)
   for j=(i+1):(n) % We choose pairs of indexes in all the possible ways.
      cont=cont + weight(i) * weight(j) * estaEntre([i j],matrizDatos); % With this subfuntion we check whether the observations are inside the band defined by observation i and j.  
      sum_weight = sum_weight + weight(i) * weight(j);               
   end
end	
   contg=cont/sum_weight;
end

function [contg]=BD3(matrizDatos, weight)
% With this function we calculate the band depth with J=3. The imput is the
% data matrix, where the rows are the observations and the columns the variables. 
[n,p]=size(matrizDatos); 
cont=zeros(1,n);      % Initialize the depth of each observation to zero.
                      % Select three observations from the sample in all the possible ways. 
sum_weight = 0;
for i=1:(n-2)
   for j=(i+1):(n-1)
      for k=(j+1):n
         cont = cont + weight(i) * weight(j) * weight(k) * estaEntre([i j k],matrizDatos);   % In this subfunction we check which observations from the sample is inside the band delimeted by observations i,j and k.           
         sum_weight = sum_weight + weight(i) * weight(j) * weight(k);
      end
   end
end	
contg=cont/sum_weight;
end

function resultado=a(v,matrizDatos) 
[n,p]=size(matrizDatos);
Z=matrizDatos;
inf=(min(Z(v,:)))';
sup=(max(Z(v,:)))';
resul=sum((and(Z'<=sup*ones(1,n),Z'>=inf*ones(1,n))));
resultado=(resul/p); % Proportion of coordinates of each observation from the sample that is inside the band delimited by pairs v of functions from the sample.
end

function [contg]=MBD(matrizDatos, weight)
% This function calculates the generalized band depth of a set of data
tic
[n,p]=size(matrizDatos); % size of the data matrix
cont=zeros(1,n); % The initial value of the depth of the observations from the matrix is zero
sum_weight = 0;
for i=1:(n-1)
    for j=(i+1):(n) % consider all possible pairs of functions
        cont = cont + weight(i) * weight(j) * a([i j],matrizDatos);  % In cont we save the generalized band depth of each observation (function)
        sum_weight = sum_weight + weight(i) * weight(j);
    end
end
contg=cont/sum_weight;
toc
end

