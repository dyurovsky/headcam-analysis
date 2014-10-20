%%% TITLE: A Scale Invariant Measure of Image Clutter
%%%
%%% USAGE: [N,K,C] = SIclutter( INFILE, EXPONENT )
%%%        - INFILE is the filename of the input image
%%%        - EXPONENT is the exponent that describes the power law
%%%          relationship between the number of segmented regions and
%%%          the scale. This exponent should be customized for one's
%%%          own specific class of images.
%%%        - N are the number of segments at each 'scale' K
%%%        - K are the scales at which segmentation is performed
%%%        - C is the clutter measure derived from N and K
%%%
%%% EXAMPLE: [N,K,C] = SIclutter( 'test.jpg', -1.3925 );
%%%
%%% NOTE: This matlab script requires the segmentation routines by
%%%       Fetzenszwalb and Huttenlocher: people.cs.uchicago.edu/~pff/segment
%%%       The executable 'segment' should be in the same directory as
%%%       this matlab script.
%%%
%%% DATE: February 10, 2011
%%%
%%% AUTHORS: Hany Farid and Mary J. Bravo
%%%
%%% CITATION:
%%%   @ARTICLE{bravo-farid08,
%%%     AUTHOR = "M.J. Bravo and H. Farid",
%%%     TITLE = "A Scale Invariant Measure of Image Clutter",
%%%     JOURNAL = "Journal of Vision",
%%%     NUMBER = "8",
%%%     VOLUME = "1",
%%%     PAGES = "1-9",
%%%     YEAR = "2008",
%%%     URL = "www.cs.dartmouth.edu/farid/publications/jov07.html"
%%%   }
%%%

function[N,K,C] = SIclutter( infile, exponent )

   K = [250 500 750 1000 1500 2000]; % segmentation parameters

   %%% LOAD, MEDIAN FILTER, AND SAVE AS .ppm
   im  = double( infile );
   %im  = double( imread( infile ) );
   im2 = zeros( size(im) );
   for z = 1 : size(im,3)
      im2(:,:,z) = medfilt2( im(:,:,z), [3 3] );
   end
   imwrite( uint8(im2), 'tmp_segin.ppm' );

   %%% SEGMENT
   sigma = 0.5;
   for k = K
      minsz = round( 0.1*k);
      cmd = sprintf('!./segment %f %d %d tmp_segin.ppm tmp_segout_%d.ppm', ...
		    sigma, k, minsz, k );
      eval( cmd );
   end

   %%% CLUTTER MEASURE
   c = 1;
   for k = K
      map  = imread( sprintf('tmp_segout_%d.ppm',k) );
      [ydim,xdim,zdim] = size( map );
      map  = reshape( map, xdim*ydim, zdim );
      map  = unique( map, 'rows' );
      N(c) = size(map,1) - 1;
      c    = c + 1;
   end
   % exponent = -1.329486; % exponent (from JOV 2007 paper)
   C = mean( N./ (K.^exponent) ); % clutter measure

   %%% REMOVE TEMP FILES
   eval( sprintf('!/bin/rm tmp_seg*.ppm') );

