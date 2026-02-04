function [grads,slews,opt] = spiralGen(fovxy,resxy,opt,plotting)
%
% This code is adapted from spiralgen_jgp_12oct.c by Jim Pipe, found at
% https://www.ismrm.org/mri_unbound/sequence.htm
%
% Changes by C. A. Baron:
%   - addition of two new spiral types:
%        sptype = 4; blipped on z to make multiple spiral planes 
%        sptype = 5; LOTUS see ISMRM 2025 abstract 1371 (Sothynathan et al) 
%
%   Inputs:
%       fovxy    in-plane field of view, in m 
%       resxy    in-plane resolution, in m 
%       opt      options structure. See code for descriptions
%

% Original comments from Pipe:
% /*********************************************
% // Spiral Generation code
% **********************************************
% // Author: Jim Pipe 
% // Date: May 2011
% // Rev: Oct 2012
% *********************************************/
% // A Subset of Relevant Literature
% //
% // Spiral Invented:
% // High-speed spiral-scan echo planar NMR imaging-I.
% // Ahn, C.B., Kim, J.H. & Cho, Z.H., IEEE Transactions on Medical Imaging, 5(1) 1986.
% //
% // Spiral Improved:
% // Fast Spiral Coronary Artery Imaging.
% // Meyer CH, Hu BS, Nishimura DG, Macovski A, Magnetic Resonance in Medicine, 28(2) 1992.
% //
% // Variable Density Spiral
% // Reduced aliasing artifacts using variable-density k-space sampling trajectories.
% // Tsai CM, Nishimura DG, Magnetic Resonance in Medicine, 43(3), 2000
% //
% // "SLOPPY" SPIRAL
% // Faster Imaging with Randomly Perturbed Undersampled Spirals and L_1 Reconstruction
% // M. Lustig, J.H. Lee, D.L. Donoho, J.M. Pauly, Proc. of the ISMRM '05
% //
% // FLORET
% // A new design and rationale for 3D orthogonally oversampled k-space trajectories
% // Pipe JG, Zwart NR, Aboussouan EA, Robison RK, Devaraj A, Johnson KO, Mag Res Med 66(5) 2011
% //
% // Distributed Spirals
% // Distributed Spirals: A New Class of 3D k-Space Trajectories
% // Turley D, Pipe JG, Magnetic Resonance in Medicine, in press (also proc of ISMRM '12)
% 
%   This function 
%   returns a single spiral arm calculated numerically
% 
%   The corresponding gradient waveforms are in gxarray and gyarray
%   spgrad_na reflects the number of gradient points to reach the end of k-space
%   spgrad_nb = spgrad_na + the number of gradient points to ramp G to zero
%   spgrad_nc = spgrad_nb + the number of gradient points to rewind k to zero
%   spgrad_nd = spgrad_nc + the number of gradient points for first moment compensation
% 
%   Assignments below indicate units of input parameters
%   All units input using kHz, msec, mT, and m!
% 
%   grad = gm exp(i theta) i.e. gm, theta are magnitude and angle of gradient
%   kloc = kr exp(i phi)   i.e. kr, phi are magnitude and angle of k-space
%   alpha = theta - phi    the angle of the gradient relative to that of k-space
%                          (alpha = Pi/2, you go in a circle
%                           alpha = 0, you go out radially)
% 
%   The variable rad_spacing determines the radial spacing
%   in units of the Nyquist distance.
%   rad_spacing = 1 gives critical sampling
%   rad_spacing > 1 gives undersampling
%   rad_spacing can vary throughout spiral generation to create variable density spirals
% 
%   KEY EQUATIONS:
%   (1) dkr/dphi = rad_spacing*Nyquist/(2 pi)
%   (2) dphi/dt = gamma gm Sin(alpha)/kr
%   (3) dkr/dt = gamma gm Cos(alpha)
% 
%   Solving (1)*(2) = (3) gives
%   (4) Tan(alpha) = (2*pi*kr)/(rad_spacing*Nyquist)
% 
% *************************************************************/
% /* Initializations */
% /************************************************************/

%%%%% Set limits
maxarray = 128000;

%%%%% Set input defaults
if nargin < 1
    fovxy = 0.220/3;
end
if nargin < 2
    resxy = 0.002;
end
if nargin < 3
    opt = [];
end
if nargin<4 || isempty(plotting)
    plotting = true;
end
if ~isfield(opt,'subrast') || isempty(opt.subrast)
    opt.subrast = 5;      %/* number of numerical cycles per gradient raster time */ 
end
if ~isfield(opt,'m_dGRast') || isempty(opt.m_dGRast)
    opt.m_dGRast = 0.01; % base raster time [msec]
end
if ~isfield(opt,'gamma') || isempty(opt.gamma)
    opt.gamma    = 42.577;    %/* typically 42.577 kHz/mT */
end
if ~isfield(opt,'gmax') || isempty(opt.gmax)
    opt.gmax     = 30;   %/* max gradient amplitude in mT/m */
end
if ~isfield(opt,'slewmax') || isempty(opt.slewmax)
    opt.slewmax  = 120; %/* max slew rate, in mT/m/msec*/
end
if ~isfield(opt,'gtype') || isempty(opt.gtype)
    % 0 = calculate through readout
    % 1 = include grad ramp-down
    % 2 = include rewinder to end at k=0
    % 3 = include first moment comp
    opt.gtype    = 1;       
end
if ~isfield(opt,'fovz') || isempty(opt.fovz)
    % Only relevant for opt.sptype == 2
    opt.fovz     = 0.256;   %    /* enter in m */ 
end
if ~isfield(opt,'resz') || isempty(opt.resz)
    opt.resz     = 0.002;   %    /* enter in m : this should be true resolution */
end
if ~isfield(opt,'arms') || isempty(opt.arms)
    opt.arms     = 1;   %    /* number of spiral interleaves*/
end
if ~isfield(opt,'sptype') || isempty(opt.sptype)
    % 0 = Archimedean
    % 1 = Cylinder DST 
    % 2 = Spherical DST
    % 3 = Fermat:Floret 
    % 4 = Archimedean with SMS CAIPI blips (CB 202301)
    % 5 = Archimedean with sinusoidal kz for SMS (CB 202301). 
      % "LOTUS": Laterally Oscillating Trajectory for Undersampling Slices
    opt.sptype    = 0;                      
end
%   /* the next 4 variables are for variable density spirals */
%   /* they create a transition in the radial spacing as the k-space radius goes from 0 to 1, i.e.*/
%   /*    0 < kr < us_0 : spacing = Nyquist distance */
%   /* us_0 < kr < us_1 : spacing increases to us_r (affected by opt.ustype)*/
%   /* us_1 < kr < 1    : spacing = us_r*/
if ~isfield(opt,'ustype') || isempty(opt.ustype)
    % rate of change in undersampling
    %     0 = linear
    %     1 = quadratic
    %     2 = hanning */
    opt.ustype   = 0; 
end
if ~isfield(opt,'us_0') || isempty(opt.us_0)
    opt.us_0    = 0; 
end
if ~isfield(opt,'us_1') || isempty(opt.us_1)
    opt.us_1    = 0; 
end
if ~isfield(opt,'us_r') || isempty(opt.us_r)
    opt.us_r    = 1; 
end
if ~isfield(opt,'slop_per') || isempty(opt.slop_per)
    %   For sloppy spirals, this lets us define periodicity in units of iteration loop time */
    %   set this to zero if you do not want sloppy spirals */
    opt.slop_per = 0; 
end
% Params for SMS   
if ~isfield(opt,'Crate') || isempty(opt.Crate)
    opt.Crate = 2; % number of simultaneous slices
end
if ~isfield(opt,'Cdz') || isempty(opt.Cdz)
    opt.Cdz = 0.1; % spacing between slices [m]
end
% Params for kz blips for CAIPI-like waveform (sptype == 4)
if ~isfield(opt,'Cf') || isempty(opt.Cf)
    opt.Cf = 0.3;  % fractional slew rate reserved for blips. 0.2 or 0.3 seems to be a good tradeoff
end
if ~isfield(opt,'CstaggerBlip') || isempty(opt.CstaggerBlip)
    opt.CstaggerBlip = 0; %// to stagger blips to try to have fewer k-space gaps
end
% Params for kz oscillations for CAIPI-like waveform (sptype == 5)
if ~isfield(opt,'CphiFact') || isempty(opt.CphiFact)
    %   period of kz wave is 1/CphiFact larger than kxy rotation period. 
    %   Good choices take many repetitions to get back to an integer to
    %   moreorless isotropically sample kz. It may also be a good idea to
    %   choose a value close to 1 so that all channels have similar freq
    %   content, which makes it easier to avoid vibrational
    %   resonances.
    opt.CphiFact = 0.618; % Equal to 1/golden ratio. 
end


%%%%% Internal variable calculation
rast     = opt.m_dGRast / opt.subrast;   %/* calculation "raster time" in msec */
if ( (opt.sptype >= 4) && (opt.Crate < 2) )
    % If only one slice, no need to account for SMS
	opt.sptype = 0;
end

%%%%% Error checking
if (opt.CstaggerBlip)
    error('Proper handling of rewinding time not tested for staggered blips')
end

%%%%% Start computations
nyquist = opt.arms/fovxy; %/* radial distance per arm to meet the Nyquist limit*/
gamrast = opt.gamma*rast; %/* gamrast*g = dk*/
dgc     = opt.slewmax*rast; %/* the most the gradients can change in 1 raster period*/
sub_gamrast = opt.subrast*gamrast;
sub_dgc     = opt.subrast*dgc;


uz=0; 
gx=0;
gy=0;
gz=0;
kx = zeros(opt.subrast*maxarray,1);
ky = zeros(opt.subrast*maxarray,1);
kz = zeros(opt.subrast*maxarray,1);
gsign = ones(opt.subrast*maxarray,1);
gxarray = zeros(maxarray,1);
gyarray = zeros(maxarray,1);
gzarray = zeros(maxarray,1);

krmax = 0.5/resxy;
kzmax = 0.5/opt.resz;
krmax2 = krmax*krmax;
kzmax2 = kzmax*kzmax;
krlim = krmax*(1.-(resxy/fovxy));

if (opt.sptype==4) 
      %/* Determine k-space spacing based on "FOV" Crate*dz. This is the
      %   typical dkz provided by a blip. 
      CkbTot = 1.0 / (opt.Crate*opt.Cdz);
      %/* Determine total k-space step provided by largest blip, which 
      %   brings you back to the starting k-space position. We design to
      %   this, then scale the smaller blips back down. */
      CkbTot = (opt.Crate - 1.0) * CkbTot;
      %// Determine duration based on basic k and grad relations [ms]
      Ct = 2*sqrt(CkbTot / (opt.gamma * opt.Cf * opt.slewmax));
      %// Find total number of points. Make total time a multiple of the twice the raster time
      Cnp = floor(Ct/(2.0*opt.m_dGRast) + 1.0) * 2 * opt.subrast;
      Ct = Cnp * rast;
      %// Find gradient change per point in array using k to grad relationship
      Cgmax = 2.0 * CkbTot / (opt.gamma * Ct);
      CGstep = Cgmax / (Cnp/2);
      %// Fill array with k-space values
      Cgb = zeros(Cnp,1);
      Cgb(1) = 0;
      for i=2:Cnp %(i=1;i<Cnp;i++) 
          if (i<=Cnp/2)
              Cgb(i) = (i-1)*CGstep;
          else
              Cgb(i) = (Cnp-i+1)*CGstep;
          end
      end
      Cbstart = -1;
      CbstartPrev = -1;
      Ccurrblip = -1;
      CcurrblipPrev = -1;
      Cadd = 0;  %// To stagger where the blips start a bit
elseif (opt.sptype == 5) 
    %// Determine k-space spacing based on "FOV" Crate*dz;  
    CkbTot = 1.0 / (opt.Crate*opt.Cdz);
    %// Determine total span to go from min to max k
    CkbTot = (opt.Crate - 1.0) * CkbTot;
    isStart = 1; 
    isStartRecord = ones(opt.subrast*maxarray,1);
    phiUnwrapped = zeros(opt.subrast*maxarray,1);
end

%/* start out spiral going radially at max slew-rate for 2 time-points */
kx(1) = 0;
ky(1) = 0;
kx(2) = gamrast*dgc;
ky(2) = 0;
kx(3) = 3*gamrast*dgc;
ky(3) = 0;

%// IF SPHERE
if (opt.sptype == 2) 
    kz(1) = kzmax;
    kz(2) = sqrt(kzmax2*(1-((kx(1)*kx(1)+ky(1)*ky(1))/krmax2))); %// stay on surface of ellipsoid
    kz(3) = sqrt(kzmax2*(1-((kx(2)*kx(2)+ky(2)*ky(2))/krmax2))); %// stay on surface of ellipsoid
end

nLoops = 0;
i = 3;
kr = kx(3);

% /******************************/
% /* LOOP UNTIL YOU HIT MAX RES */
% /******************************/
while ((kr <= krlim) && (i < opt.subrast*maxarray-1) ) 
    if (nLoops > 10*opt.subrast*maxarray)
        error('Spiral gen failure')
    end
    if (i<2)
        error('Spiral gen failure on rewinding time')
    end

    % /**************************/
    % /*** STEP 1:  Determine the direction (ux,uy) of the gradient at ~(i+0.5) */
    % /**************************/
    %    /* calculate dk/rast = opt.gamma G*/
    kmx = 1.5*kx(i) - 0.5*kx(i-1);
    kmy = 1.5*ky(i) - 0.5*ky(i-1);
    kmr = sqrt(kmx*kmx + kmy*kmy);

    % /////////////////////////////
    % // Start rad_spacing logic //
    % /////////////////////////////
    rnorm = 2*resxy*kmr; %/* the k-space radius, normalized to go from 0 to 1 */
    
    %/* determine the undersample factor */
    if (rnorm <= opt.us_0)
        rad_spacing = 1;
    elseif (rnorm < opt.us_1) 
        us_i = (rnorm-opt.us_0)/(opt.us_1 - opt.us_0); %/* goes from 0 to 1 as rnorm goes from us_0 to us_1*/
        if (opt.ustype == 0) 
            %/* linearly changing undersampling*/
            rad_spacing = 1. + (opt.us_r - 1.)*us_i;
        elseif (opt.ustype == 1) 
            %/* quadratically changing undersampling*/
            rad_spacing = 1. + (opt.us_r - 1.)*us_i*us_i;
        elseif (opt.ustype == 2) 
            %/* Hanning-type change in undersampling */
            rad_spacing = 1. + (opt.us_r - 1.)*0.5*(1.-cos(us_i*M_PI));
        end
    else 
        rad_spacing = opt.us_r;
    end

    %/* Undersample spiral for Spherical-Distributed Spiral */
    if (opt.sptype == 2) 
        if (rnorm < 1.0)
            rad_spacing = min(opt.fovz/opt.resz, rad_spacing/sqrt(1.0 - (rnorm*rnorm)));
        else
            rad_spacing = opt.fovz/opt.resz;
        end
    end

    %/* MAKE FERMAT SPIRAL FOR FLORET*/
    if (opt.sptype == 3 && rnorm > 0) 
        rad_spacing = rad_spacing / rnorm;
    end

    %/* Sloppy Spirals - add variability to rad_spacing for reduced aliasing coherence */ 
    % // A couple different options here are commented out
    % // Lots of ways to be sloppy
    if (opt.slop_per > 0) 
        % //      rad_spacing = MAX(1., (rad_spacing + ((rad_spacing-1.)*sin(2.*M_PI*(double)(i)/opt.slop_per))));
        % //      rad_spacing += (rad_spacing-1.)*sin(2.*M_PI*opt.slop_per*atan2(ky(i),kx(i)));
        rad_spacing = rad_spacing + (rad_spacing-1)*sin(2*pi*opt.slop_per*rnorm);
    end

    % ///////////////////////////
    % // End rad_spacing logic //
    % ///////////////////////////

    %/* See the Key Equation 4 at the beginning of the code */
    alpha = atan(2*pi*kmr/(rad_spacing*nyquist));
    phi = atan2(kmy,kmx);
    theta = phi + alpha;

    ux = cos(theta);
    uy = sin(theta);

    % // IF SPHERICAL DST
    % // u dot km is zero if moving on a sphere (km is radial, u is tangential,
    % // thus km stays on the sphere)
    % // We are on an ellipsoid, but can normalize u and km by krmax and kzmax to make this work
    % // The final gradient vector (ux uy uz) will be tangential to the sphere
    if (opt.sptype == 2) 
        kmz = 1.5*kz(i) - 0.5*kz(i-1);
        uz = -((ux*kmx + uy*kmy)/krmax2)*(kzmax2/kmz);
        umag = sqrt(ux*ux + uy*uy + uz*uz);
        ux = ux/umag;
        uy = uy/umag;
        uz = uz/umag;
        gz = (kz(i) - kz(i-1))/gamrast;
    elseif (opt.sptype == 4)
        % Set gradient for pre-defined CAIPI blip
        gz = (kz(i) - kz(i-1))/gamrast;
        %// Start a blip whenever kx goes from pos to neg
        if ( (kx(i-1) > 0) && (kx(i) < 0) )
            if ( (Cbstart > 0) && (i >= Cbstart) && (i < Cbstart + Cnp ) )
                error('CAIPI blips overlapping')
            end
            CbstartPrev = Cbstart;
		    Cbstart = i + Cadd;
            CcurrblipPrev = Ccurrblip;
            Ccurrblip = Ccurrblip + 1;
            if (Ccurrblip > opt.Crate)
                Ccurrblip = 1;
            end
            if opt.CstaggerBlip
                Cadd = Cadd + round(Cnp/opt.Crate);
                if Cadd > Cnp
                    Cadd = 0;
                end
            end
        end
        %// Check if in blip
        if ( (Cbstart > 0) && (i >= Cbstart) && (i < Cbstart + Cnp ) ) 
            %// Scale blip and set polarity 
            if (Ccurrblip == 0)
                %// First blip is scaled differently to make kz symmetric
                Cfact = -0.5;
            elseif (Ccurrblip < opt.Crate)
                Cfact = 1/(opt.Crate-1);
            else
                Cfact = -1;
            end
            gznext = Cfact*Cgb(i-Cbstart+1); 
        else
            gznext = 0.0;
        end
    elseif (opt.sptype == 5)
        %/* Find unwrapped dphi. */
        dphi = atan2(ky(i),kx(i)) - atan2(ky(i-1),kx(i-1));
        if abs(dphi)>3*pi/4
            if (dphi<0)
                dphi = dphi + 2*pi;
            else
                dphi = dphi - 2*pi;
            end
        end
        phiUnwrapped(i) = phiUnwrapped(i-1) + dphi;
        % Our target is kz = CkbTot/2*cos(opt.CphiFact*phi). 
        %   Thus, dk/dphi = CkbTot/2*opt.CphiFact*sin(opt.CphiFact*phi)
        %   From comments at top of file, dphi/dt = gamma*Gxy*sin(alpha)/kmr
        %   Multiplying these equations: dk/dt = gamma*Gxy*sin(alpha)/kmr*CkbTot/2*opt.CphiFact*sin(opt.CphiFact*phi)
        %   Recognizing that gamma*Gz = dk/dt, using ux,uy,uz for Gx,Gy,Gz, and rearranging yields:
        uz = sqrt(ux^2+uy^2)*opt.CphiFact*CkbTot/2/kmr*...
            sin(opt.CphiFact*phiUnwrapped(i))*sin(alpha);
        % If we just use the above uz, the kz span will not be centered.
        % So, we use the first half-period of cos(opt.CphiFact*phi) as a
        % "prephasor" to get to the edge of the desired span. We can do
        % this by scaling uz by 0.5 during the first pi radians of opt.CphiFact*phi
        isStartRecord(i) = isStart;
        if (isStart)
            uz = 0.5*uz;
            if (opt.CphiFact*phiUnwrapped(i) >= pi)
                isStart = 0;
            end
        end
        %// Normalize the unit vector
        umag = sqrt(ux*ux + uy*uy + uz*uz);
        ux = ux/umag;
        uy = uy/umag;
        uz = uz/umag;
        gz = (kz(i) - kz(i-1))/gamrast;
    end

    % /**************************/
    % /*** STEP 2: Find largest gradient magnitude with available slew */
    % /**************************/

    %/* Current gradient*/
    gx = (kx(i) - kx(i-1))/gamrast;
    gy = (ky(i) - ky(i-1))/gamrast;

    % /*
    % // solve for gm using the quadratic equation |gm u - g| = dgc
    % // which is
    % //   (gm u - g)(gm u* - g*) = dgc^2
    % // which gives
    % //   gm^2 (u u*) - gm (g u* + u g*) + g g* - dgc^2 = 0
    % 
    % // Replacing u u* with 1 (i.e. u is a unit vector) and
    % // replacing (g u* + u g*) with 2 Real[g u*]
    % // this is
    % //   gm^2 + gm (2 b) + c = 0
    % // giving
    % //   gm = -b +/- Sqrt(b^2 - c)
    % // The variable "term" = (b^2 - c) will be positive if we can meet the desired new gradient 
    % */
    if (opt.sptype == 4)
        %// Only allow slew not reserved for blips. Ignore gz
        %// keep slew reduced for whole waveform. Could be more efficient, but this is easy
        term = dgc*dgc*(1-opt.Cf^2) - (gx*gx + gy*gy) + (ux*gx + uy*gy)*(ux*gx + uy*gy);
    else
        term = dgc*dgc - (gx*gx + gy*gy + gz*gz) + (ux*gx + uy*gy + uz*gz)*(ux*gx + uy*gy + uz*gz);
    end

    if (term >= 0) 
        % // Slew constraint is met! Now assign next gradient and then next k value
        % // NOTE gsign is +1 or -1
        % //   if gsign is positive, we are using slew to speed up (increase gm) as much as possible
        % //   if gsign is negative, we are using slew to slow down (decrease gm) as much as possible
        if (opt.sptype == 4)
            %// Account for fixed blip gradient contributing to net grad
            %// keep xy max grad reduced for whole waveform. Could be more efficient, but this is easy
            gm  = min((ux*gx + uy*gy) + gsign(i)*sqrt(term),sqrt(opt.gmax^2-Cgmax^2));
        else
            gm  = min((ux*gx + uy*gy + uz*gz) + gsign(i)*sqrt(term),opt.gmax);
        end
        gx = gm*ux;
        gy = gm*uy;

        kx(i+1) = kx(i) + gx*gamrast;
        ky(i+1) = ky(i) + gy*gamrast;

        %// If SPHERE
        if (opt.sptype == 2)
            kz(i+1) = sqrt(kzmax2*(1.-((kx(i+1)*kx(i+1)+ky(i+1)*ky(i+1))/krmax2))); %// stay on surface of ellipsoid
        elseif (opt.sptype == 4)
            kz(i+1) = kz(i) + gznext*gamrast;
        elseif (opt.sptype == 5)
            gz = gm*uz;
            kz(i+1) = kz(i) + gz*gamrast;
        end

        i = i+1;
    else 
        % // We can't go further without violating the slew rate
        % // This means that we've sped up too fast to turn here at the desired curvature
        % // We are going to iteratively go back in time and slow down, rather than speed up, at max slew
        % // Here we'll keep looking back until gsign is positive, then add another negative gsign, just far enough to make the current corner
      while ((i>4) && (gsign(i-1) == -1)) 
          i = i-1;
      end
      gsign(i-1) = -1;
      i = i-2;
      if (opt.sptype == 4) && (i<=Cbstart)
          % We rewound past a blip start
          Cbstart = CbstartPrev;
          Ccurrblip = CcurrblipPrev;
      elseif (opt.sptype == 5)
          isStart = isStartRecord(i);
      end
    end

    kr = sqrt(kx(i)*kx(i) + ky(i)*ky(i));
    
    nLoops = nLoops + 1;

end % End main kr while loop

i_end = i;

% //********************************************
% // DONE LOOPING FOR SAMPLING PORTION
% // recast k to g while subsampling by opt.subrast  
% //********************************************
% TODO: vectorize this
gxsum = 0;
gysum = 0;
gzsum = 0;
for j = 1:floor(i_end/opt.subrast) %(j=1;j<=(i_end/opt.subrast);j++) 
    i1 = j*opt.subrast + 1;
    i0 = (j-1)*opt.subrast + 1;
    gxarray(j) = ( (kx(i1)-kx(i0))/sub_gamrast );
    gyarray(j) = ( (ky(i1)-ky(i0))/sub_gamrast );
    gzarray(j) = ( (kz(i1)-kz(i0))/sub_gamrast );
    gxsum = gxsum + gxarray(j);
    gysum = gysum + gyarray(j);
    gzsum = gzsum + gzarray(j);
end
spgrad_na = j;

%// recalculate these ending gradient points
gm = sqrt(gxarray(spgrad_na-1)*gxarray(spgrad_na-1) +...
        gyarray(spgrad_na-1)*gyarray(spgrad_na-1) +...
        gzarray(spgrad_na-1)*gzarray(spgrad_na-1));
ux = gxarray(spgrad_na-1)/gm;
uy = gyarray(spgrad_na-1)/gm;
uz = gzarray(spgrad_na-1)/gm;

% //**************************************************
% // NOW, if requested via gtype, go to g=0 and k=0
% // I've tried other ways to be faster, can't find them
% //**************************************************
  
% // first we'll ramp gradients to zero
% // note {ux,uy} is still pointing in the gradient direction
% TODO: vectorize
if (opt.gtype > 0) 
    gz_sum_ramp = 0;
    while ((gm > 0) && (j < maxarray)) 
        gm = max(0,gm - sub_dgc);
        gxarray(j) = gm*ux;
        gyarray(j) = gm*uy;
        gzarray(j) = gm*uz;
        gxsum = gxsum + gxarray(j);
        gysum = gysum + gyarray(j);
        gzsum = gzsum + gzarray(j);
        gz_sum_ramp = gz_sum_ramp + gzarray(j);
        j = j+1;
    end
end
spgrad_nb = j;

% // now point gradient towards the k-space origin
% // {ux,uy} will be a unit vector in that direction
if (opt.gtype > 1) 
    %     /* NOTE: spherical needs a prephaser not a rewinder 
    %      * so just rewind x and y in that case */
    gsum = sqrt(gxsum*gxsum + gysum*gysum + gzsum*gzsum);
    if (opt.sptype == 2 ) 
        gsum = sqrt(gxsum*gxsum + gysum*gysum + gz_sum_ramp*gz_sum_ramp);
    end
    gsum0 = gsum;
    ux = -gxsum/gsum;
    uy = -gysum/gsum;
    uz = -gzsum/gsum;
    if (opt.sptype == 2) 
        uz = -gz_sum_ramp/gsum;
    end
    gsum_ramp = 0.5*gm*(gm/sub_dgc); %/* this is *roughly* how much the area changes if we ramp down the gradient NOW*/
                                     %/* this value is zero right now (gm = 0), but it will make sense below */
    %// increase gm while we can
    while ((gsum_ramp < gsum) && (j < maxarray)) 
        gm = min(opt.gmax,gm+sub_dgc);
        gxarray(j) = gm*ux;
        gyarray(j) = gm*uy;
        gzarray(j) = gm*uz;
        gsum = gsum - gm;
        j = j+1;
        gsum_ramp = 0.5*gm*(gm/sub_dgc); %/* see - now this makes sense; this tells us when to start ramping down */
    end

    % // We've overshot it by a tiny bit, but we'll fix that later
    % // Ramp down for now
    while ((gm > 0) && (j < maxarray)) 
      gm = max(0,gm-sub_dgc);
      gxarray(j) = gm*ux;
      gyarray(j) = gm*uy;
      gzarray(j) = gm*uz;
      gsum = gsum - gm;
      j = j+1;
    end
    spgrad_nc = j;

    %// OK - gm is zero, but gsum is probably not EXACTLY zero. Now scale the rewinder to make the sum exactly zero
    gradtweak = gsum0/(gsum0-gsum);
    for j = spgrad_nb+1:spgrad_nc %(j=(spgrad_nb); j<(spgrad_nc); j++) 
        gxarray(j) = (gradtweak)*gxarray(j);
        gyarray(j) = (gradtweak)*gyarray(j);
        gzarray(j) = (gradtweak)*gzarray(j);
    end
end

gxarray = gxarray(1:j);
gyarray = gyarray(1:j);
gzarray = gzarray(1:j);

grads = cat(2, gxarray, gyarray, gzarray);
slews = diff(grads,1,1)/opt.m_dGRast;


%% Plotting
if plotting
    figure;
    n1 = 2;
    n2 = 4;
    subplot(n1,n2,1); 
    plot(grads);
    title('gradients')
    
    subplot(n1,n2,2)
    kx = cumsum(gxarray)*opt.gamma*opt.m_dGRast;
    ky = cumsum(gyarray)*opt.gamma*opt.m_dGRast;
    kz = cumsum(gzarray*opt.gamma*opt.m_dGRast);
    plot(kx, ky)
    if opt.sptype == 5
        zc = [abs(diff(sign(kz)))>0;false];
        kx_a = kx(zc);
        ky_a = ky(zc);
        hold('all');
        plot(kx_a,ky_a,'o')
    end
    title('k-space traj xy with kz zero crossings')
    
    subplot(n1,n2,3) 
    plot3(cumsum(gxarray*opt.gamma*opt.m_dGRast), cumsum(gyarray*opt.gamma*opt.m_dGRast), cumsum(gzarray*opt.gamma*opt.m_dGRast))
    title('k-space traj xyz')
    
    subplot(n1,n2,4) 
    % plot(cumsum(gxarray*opt.gamma*opt.m_dGRast))
    % hold('all')
    % plot(cumsum(gyarray*opt.gamma*opt.m_dGRast))
    plot(cumsum(gzarray*opt.gamma*opt.m_dGRast))
    title('k-space traj z')

    subplot(n1,n2,n2+1); 
    plot(sqrt(sum(slews.^2,2)));
    title('net slew')
    
    subplot(n1,n2,n2+2); 
    plot(sqrt(sum(grads.^2,2)));
    title('net grad')
    
    % Freq analysis notes:
    %   Gmax and slew has the largest effect. 
    %   Undersampling and variable density has little effect
    %   CAIPI grads are negligible compared to others.
    subplot(n1,n2,n2+3)
    fmaxp = 3000;
    Nf = 10*length(gxarray);
    fmax = 1000*0.5/opt.m_dGRast;  % Hz
    psd = fft([gxarray,gyarray,gzarray],Nf,1);
    psd = psd(1:ceil(Nf/2),:).*conj(psd(1:ceil(Nf/2),:));
    f = linspace(0,fmax,size(psd,1))';
    plot(f,psd)
    xlim([0,fmaxp])
    cent = sum(f.*psd,1)./sum(psd,1);
    [~, pk] = max(psd(:,1)); pk = f(pk);
    bw = find(psd(:,1) > max(psd(:,1))/2);
    bw = f(bw(end)) - f(bw(1));
    text(0.5*fmaxp, 0.9*max(psd(:,1)),... 
        sprintf('cent = %d Hz\npeak = %d Hz\nbw = %d Hz\nrate %d\nvd %.2f',round(cent(1)),round(pk(1)),round(bw),opt.us_r,opt.us_1))
end


