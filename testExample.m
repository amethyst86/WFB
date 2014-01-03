% testing weighted functional boxplot using synthetic data 

close all
clear all

% parameter settings for weighted functional boxplots 
wfbParams.factor = 1.5;  
wfbParams.barColor = 'b';
wfbParams.centerColor = 'm';
wfbParams.showOut = true;
wfbParams.method = 'MBD';

% set sigma for gaussian function
gaussianSigma = 30;

% generate random ages for each curve
ages = min( max( round(180 * rand( 30, 1 )), 1 ), 180 ); 
ages = sort(ages);
nBins = ceil(max(ages)/10)*10;
cmap = jet(nBins);
% age changing range
ageRange = [0 200];

% generate curves
% each colomn of the yAxis matrix is corresponding to a curve
xAxis = (0:0.01:1)';
yAxis = zeros( length(xAxis), length(ages) );

% plot the generated curves
figure, hold on
for iI = 1:length(ages)
    yAxis( :, iI ) = 500 * ( 1 + sin( xAxis*3.14*2 + 1.57*iI/5 ) ) + ages(iI)*2;
    plot( xAxis, yAxis( :, iI ), 'Color', cmap(ages(iI), :), 'LineWidth', 2 );
end
hold off
caxis( [0 nBins] );
h = colorbar( 'peer', gca );
set( get( h, 'ylabel' ), 'String', 'Age');

% built atlas at the age of each curve
estimatedAge = ages;
maxAge = max(ages);
figure;
for id = 1:length(estimatedAge)
  clf;
  ageId = estimatedAge( id );
  
  % using boundary reflection to mitigate the boundary bias
  weightUnnormalized = gaussianWeighting( ages, ageRange(1), ageRange(2), gaussianSigma, ageId );
  % normalize the weights
  weights = weightUnnormalized / sum( weightUnnormalized ) ;
  
  % plot the weighted functional boxplot
  subplot(2, 1, 1), [depthTmp, medianCurveId] = wfbplot( yAxis, xAxis, weights, wfbParams );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename );

  % plot the age histogram
  subplot(2, 1, 2), hist(ages, ageRange(2)-ageRange(1)+1); xlim([ageRange(1), ageRange(2)]);
  hold on
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.25, textString, 'Color', 'm', 'FontSize', 15, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d (months)', ageId );
  text( ageId-2, -0.35, textString, 'Color', 'r', 'FontSize', 15, 'FontWeight', 'Bold');
  [sort_age, sort_Id] = sort(ages);
  plot(sort_age, weightUnnormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
   
  drawnow
end
