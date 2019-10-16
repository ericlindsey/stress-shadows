function [E] = computeStrainMeade07(sx,sy,sz,x,y,z,pr,ss,ts,ds)
% TRIANGLESTRAINS calculates strains due to slip on a triangular dislocation
% in an elastic half space utilizing the symbolically differentiated
% displacement gradient tensor derived from the expressions for
% the displacements due to an angular dislocation in an elastic half
% space (Comninou and Dunders, 1975).
%
% INPUT:
%  sx - x-coordinates of observation points (east)
%  sy - y-coordinates of observation points (north)
%  sz - z-coordinates of observation points (positive down)
%  x  - x-coordinates of triangle vertices.
%  y  - y-coordinates of triangle vertices.
%  z  - z-coordinates of triangle vertices.
%  pr - Poisson's ratio
%  ss - strike slip displacement (positive for right-lateral slip)
%  ts - tensile slip displacement
%  ds - dip slip displacement (thrust for dip<90)
%
% OUTPUT:
%
%  E  - structure containing the strains
%       E.xx, E.yy, E.zz, E.xy, E.xz, E.yz)
%       (x is for east, y north, and z up)
%
% SEE ALSO: unicycle

% This paper should and related code should be cited as:
% Brendan J. Meade, Algorithms for the calculation of exact
% displacements, strains, and stresses for Triangular Dislocation
% Elements in a uniform elastic half space, Computers &
% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
%
% Use at your own risk and please let me know of any bugs/errors.
%
% Copyright (c) 2006 Brendan Meade
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% HISTORY: Jan 26, 2012, Sylvain Barbot. Factorize repetitive terms in the
% strain calculation. Modifications result in acceleration by a factor of 19.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

% Calculate the slip vector in XYZ coordinates
normVec=cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)],[x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec=normVec./norm(normVec);

% choose upward-pointing normal vectors
if (normVec(3) < 0)
    normVec=-normVec;
    [x(2) x(3)]=swap(x(2), x(3));
    [y(2) y(3)]=swap(y(2), y(3));
    [z(2) z(3)]=swap(z(2), z(3));
end
strikeVec=[-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec=cross(normVec, strikeVec);
slipComp=[ss ds ts];
slipVec=[strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);

% Solution vectors
E.xx=zeros(size(sx));
E.yy=zeros(size(sx));
E.zz=zeros(size(sx));
E.xy=zeros(size(sx));
E.xz=zeros(size(sx));
E.yz=zeros(size(sx));

% Add a copy of the first vertex to the vertex list for indexing
x(4)=x(1);
y(4)=y(1);
z(4)=z(1);

for iTri = 1:3
    % Calculate strike and dip of current leg
    strike=180.0/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
    segMapLength=sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
    [rx ry]=RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
    dip=180.0/pi*(atan2(z(iTri+1)-z(iTri), rx));
    
    if dip>=0
        beta=pi/180.0*(90.0-dip);
        if pi/2<beta
            beta=pi/2-beta;
        end
    else
        beta=-pi/180.0*(90.0+dip);
        if -pi/2>beta
            beta=pi/2-abs(beta);
        end
    end
    ssVec=[cos(strike/180.0*pi) sin(strike/180.0*pi) 0];
    tsVec=[-sin(strike/180.0*pi) cos(strike/180.0*pi) 0];
    dsVec=cross(ssVec, tsVec);
    lss=dot(slipVec, ssVec);
    lts=dot(slipVec, tsVec);
    lds=dot(slipVec, dsVec);
    
    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
        % First angular dislocation
        [sx1 sy1]=RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        [a11 a22 a33 a12 a13 a23]=advs(sx1,sy1,sz-z(iTri),z(iTri),beta,pr,lss,lts,lds);
        
        % Second angular dislocation
        [sx2 sy2]=RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        [b11 b22 b33 b12 b13 b23]=advs(sx2,sy2,sz-z(iTri+1),z(iTri+1),beta,pr,lss,lts,lds);
        
        % Rotate tensors to correct for strike
        bxx=a11-b11;
        byy=a22-b22;
        bzz=a33-b33;
        bxy=a12-b12;
        bxz=a13-b13;
        byz=a23-b23;
        
        g=pi/180*strike;
        e11n=(cos(g)*bxx-sin(g)*bxy)*cos(g)-(cos(g)*bxy-sin(g)*byy)*sin(g);
        e12n=(cos(g)*bxx-sin(g)*bxy)*sin(g)+(cos(g)*bxy-sin(g)*byy)*cos(g);
        e13n=cos(g)*bxz-sin(g)*byz;
        e22n=(sin(g)*bxx+cos(g)*bxy)*sin(g)+(sin(g)*bxy+cos(g)*byy)*cos(g);
        e23n=sin(g)*bxz+cos(g)*byz;
        e33n=bzz;
        
        % Add the strains from current leg
        E.xx=E.xx+e11n;
        E.yy=E.yy+e22n;
        E.zz=E.zz+e33n;
        E.xy=E.xy+e12n;
        E.xz=E.xz-e13n;
        E.yz=E.yz-e23n;
    end
end
end

function [a b] = swap(a, b)
    % Swap two values
    temp=a;
    a=b;
    b=temp;
end

function [xp yp] = RotateXyVec(x, y, alpha)
    % Rotate a vector by an angle alpha
    x=x(:);
    y=y(:);
    alpha=pi/180*alpha;
    xp=cos(alpha).*x-sin(alpha).*y;
    yp=sin(alpha).*x+cos(alpha).*y;
end

function [e11 e22 e33 e12 e13 e23] = advs(y1,y2,y3,a,b,nu,B1,B2,B3)
    % These are the strains in a uniform elastic half space due to slip
    % on an angular dislocation.  They were calculated by symbolically
    % differentiating the expressions for the displacements (Comninou and
    % Dunders, 1975, with typos noted by Thomas 1993) then combining the
    % elements of the displacement gradient tensor to form the strain tensor.
    
    t1=1-nu;
    t2=2*t1;
    t3=y1.^2;
    t4=1./t3;
    t6=(y2.^2);
    t9=1./(1+t6.*t4);
    t10=y2.*t4.*t9;
    t12=cos(b);
    t13=y1.*t12;
    t14=sin(b);
    t15=y3.*t14;
    t16=t13-t15;
    t17=t16.^2;
    t18=0.1e1./t17;
    t19=y2.*t18;
    t22=0.1e1./(0.1e1+t6.*t18);
    t24=t19.*t12.*t22;
    t25=y3.^2;
    t26=t3+t6+t25;
    t27=sqrt(t26);
    t28=0.1e1./t27;
    t29=y2.*t28;
    t31=t6.*t12;
    t32=y1.*t16+t31;
    t33=0.1e1./t32;
    t34=t14.*t33;
    t37=y2.*t27;
    t38=t32.^2;
    t39=0.1e1./t38;
    t40=t14.*t39;
    t41=0.2e1.*t13;
    t47=t14.^2;
    t48=t47.*t39;
    t51=0.1e1./(0.1e1+t6.*t26.*t48);
    t52=(t29.*t34.*y1-t37.*t40.*(t41-t15)).*t51;
    t53=2.*a;
    t54=y3+t53;
    t55=t54.*t14;
    t56=t13+t55;
    t57=t56.^2;
    t58=0.1e1./t57;
    t59=y2.*t58;
    t62=0.1e1./(0.1e1+t6.*t58);
    t64=t59.*t12.*t62;
    t65=t54.^2;
    t66=t3+t6+t65;
    t67=sqrt(t66);
    t68=0.1e1./t67;
    t69=y2.*t68;
    t71=y1.*t56+t31;
    t72=0.1e1./t71;
    t73=t14.*t72;
    t76=y2.*t67;
    t77=t71.^2;
    t78=0.1e1./t77;
    t79=t14.*t78;
    t85=t47.*t78;
    t88=0.1e1./(0.1e1+t6.*t66.*t85);
    t89=(t69.*t73.*y1-t76.*t79.*(t41+t55)).*t88;
    t91=t2.*((2.*t10)-t24+t52-t64+t89);
    t92=t27-y3;
    t93=0.1e1./t92;
    t94=t28.*t93;
    t95=t67+y3+t53;
    t96=0.1e1./t95;
    t97=t68.*t96;
    t98=t94+t97;
    t99=y2.*t98;
    t100=y1.*y2;
    t102=0.1e1./t27./t26;
    t103=t102.*t93;
    t104=t103.*y1;
    t105=0.1e1./t26;
    t106=t92.^2;
    t107=0.1e1./t106;
    t108=t105.*t107;
    t109=t108.*y1;
    t111=0.1e1./t67./t66;
    t112=t111.*t96;
    t113=t112.*y1;
    t114=0.1e1./t66;
    t115=t95.^2;
    t116=0.1e1./t115;
    t117=t114.*t116;
    t118=t117.*y1;
    t119=-t104-t109-t113-t118;
    t120=t100.*t119;
    t121=y2.*t12;
    t122=t28.*t14;
    t124=t122.*y1-0.1e1;
    t126=y1.*t14;
    t128=t27-t126-y3.*t12;
    t129=0.1e1./t128;
    t131=t27.*t14;
    t132=t131-y1;
    t133=t132.*t102;
    t134=t129.*y1;
    t136=t132.*t28;
    t137=t128.^2;
    t138=0.1e1./t137;
    t139=t28.*y1;
    t140=t139-t14;
    t141=t138.*t140;
    t143=t68.*t14;
    t145=t143.*y1-0.1e1;
    t146=t145.*t68;
    t147=t54.*t12;
    t148=t67-t126+t147;
    t149=0.1e1./t148;
    t151=t67.*t14;
    t152=t151-y1;
    t153=t152.*t111;
    t154=t149.*y1;
    t156=t152.*t68;
    t157=t148.^2;
    t158=0.1e1./t157;
    t159=t68.*y1;
    t160=t159-t14;
    t161=t158.*t160;
    t163=t124.*t28.*t129-t133.*t134-t136.*t141+t146.*t149-t153.*t154-t156.*t161;
    t166=0.1e1./pi;
    t168=1./t1;
    t171=-t2;
    t172=2.*nu;
    t173=1-t172;
    t174=t171.*t173;
    t175=t10-t64+t89;
    t176=cot(b);
    t177=t176.^2;
    t178=t175.*t177;
    t180=t173.*y2;
    t181=t180.*t116;
    t182=a.*t68;
    t183=0.1e1-t172-t182;
    t184=t183.*t176;
    t185=y1.*t96;
    t186=nu+t182;
    t187=t185.*t186;
    t188=t184-t187;
    t192=a.*t111;
    t193=y1.*t176;
    t194=t192.*t193;
    t196=t3.*t116;
    t197=t186.*t68;
    t199=t3.*t96;
    t200=t199.*t192;
    t201=t194-t96.*t186+t196.*t197+t200;
    t204=t180.*t12;
    t205=t176.*t158;
    t206=t12+t182;
    t207=t206.*t160;
    t210=t12.*t176;
    t211=t180.*t210;
    t212=t149.*a;
    t213=t111.*y1;
    t214=t212.*t213;
    t217=y3+a;
    t218=a.*y2.*t217;
    t219=t66.^2;
    t221=0.1e1./t67./t219;
    t222=t176.*t221;
    t223=t222.*y1;
    t225=0.3e1.*t218.*t223;
    t226=y2.*t217;
    t227=t226.*t111;
    t228=-t173;
    t229=t228.*t176;
    t230=t172+t182;
    t231=t185.*t230;
    t232=a.*y1;
    t233=t232.*t114;
    t234=t229+t231+t233;
    t235=t96.*t234;
    t238=t226.*t114;
    t239=t116.*t234;
    t242=t96.*t230;
    t243=t230.*t68;
    t245=a.*t114;
    t246=a.*t3;
    t247=0.1e1./t219;
    t249=0.2e1.*t246.*t247;
    t253=t12.*t149;
    t255=t67.*t12+y3+t53;
    t256=t173.*t12;
    t257=t256-t182;
    t262=t255.*t257.*t176+t2.*t152.*t12;
    t264=a.*t54;
    t267=t253.*t262-t264.*t210.*t114;
    t268=t149.*t267;
    t271=t226.*t68;
    t272=t158.*t267;
    t275=t68.*t149;
    t276=t12.*t158;
    t279=t68.*t12;
    t283=t255.*a;
    t290=t264.*t12;
    t291=t176.*t247;
    t298=t174.*t178-t181.*t188.*t68.*y1+t180.*t96.*t201-t204.*t205.*t207-t211.*t214-t225-t227.*t235.*y1-t238.*t239.*y1+t226.*t97.*(t242-t196.*t243-t200+t245-t249)-t227.*t268.*y1-t271.*t272.*t160+t226.*t275.*(-t276.*t262.*t160+t253.*(t279.*y1.*t257.*t176+t283.*t213.*t176+t2.*t145.*t12)+0.2e1.*t290.*t291.*y1);
    t305=t159.*t96;
    t306=t140.*t129;
    t307=t160.*t149;
    t308=t306+t307;
    t310=t139.*t93+t305-t12.*t308;
    t312=y1.*t98;
    t316=t28.*t129;
    t320=t16.*t132;
    t321=t102.*t129;
    t322=t321.*y1;
    t324=t28.*t138;
    t325=t324.*t140;
    t331=t56.*t152;
    t332=t111.*t149;
    t333=t332.*y1;
    t335=t68.*t158;
    t336=t335.*t160;
    t338=t228.*t310+0.2e1.*t312+t3.*t119+t12.*t132.*t316+t16.*t124.*t316-t320.*t322-t320.*t325+t12.*t152.*t275+t56.*t145.*t275-t331.*t333-t331.*t336;
    t342=t2.*t177;
    t343=t342+nu;
    t344=t343.*t68;
    t347=(t342+0.1e1).*t12;
    t351=t173.*t116;
    t354=nu.*t54;
    t355=t176.*t68;
    t358=t228.*y1.*t176+t354-a+t232.*t355+t199.*t186;
    t359=t358.*t68;
    t362=t173.*t96;
    t365=t176.*t111;
    t366=t246.*t365;
    t368=t3.*y1;
    t375=t173.*t176;
    t377=a.*t152;
    t378=0.1e1./t12;
    t379=t68.*t378;
    t381=t56.*t12-t377.*t379;
    t382=t158.*t381;
    t385=t12.^2;
    t388=t111.*t378;
    t389=t388.*y1;
    t395=a.*t217.*t365;
    t396=t217.*t176;
    t397=t396.*t221;
    t399=0.3e1.*t246.*t397;
    t400=t217.*t116;
    t402=t173.*y1.*t176;
    t403=t402+a;
    t405=t3.*t68;
    t408=t172+t68.*t403-t405.*t242-t246.*t111;
    t409=t408.*t68;
    t412=t217.*t96;
    t413=t111.*t403;
    t415=t68.*t173;
    t416=t415.*t176;
    t422=t116.*t230;
    t425=t96.*a;
    t427=t192.*y1;
    t434=t12.*t14;
    t438=t2.*t12;
    t439=t255.*t149;
    t441=0.1e1+t182.*t378;
    t442=t439.*t441;
    t443=t438-t442;
    t445=-t434+t232.*t54.*t111.*t378+t156.*t443;
    t446=t158.*t445;
    t458=t279.*t154.*t441;
    t459=t255.*t158;
    t461=t459.*t441.*t160;
    t462=t439.*a;
    t463=t462.*t389;
    t469=t173.*(t344.*t185-t347.*t307)-t351.*t359.*y1+t362.*(t229+a.*t176.*t68-t366+0.2e1.*t187-t368.*t116.*t197-t368.*t96.*t192)+t375.*t382.*t160-t375.*t149.*(t385-a.*t145.*t379+t377.*t389)-t395+t399-t400.*t409.*y1+t412.*(-t413.*y1+t416-0.2e1.*t159.*t242+t368.*t111.*t242+t368.*t114.*t422+t368.*t247.*t425-0.2e1.*t427+0.3e1.*a.*t368.*t221)-t396.*t446.*t160+t396.*t149.*(t264.*t388-0.3e1.*t246.*t54.*t221.*t378+t146.*t443-t153.*t443.*y1+t156.*(-t458+t461+t463));
    t475=y2.*t14;
    t480=y2.*t116;
    t481=0.1e1+t182;
    t482=t481.*t68;
    t484=t480.*t482.*y1;
    t485=y2.*t96;
    t486=t485.*t427;
    t487=t158.*t206;
    t488=t487.*t160;
    t489=t121.*t488;
    t490=t121.*t149;
    t491=t490.*t427;
    t494=t245+t96;
    t495=t111.*t494;
    t497=t226.*t495.*y1;
    t498=a.*t247;
    t499=t498.*y1;
    t501=t116.*t68;
    t502=t501.*y1;
    t504=t68.*(-0.2e1.*t499-t502);
    t506=t226.*t12;
    t507=t439.*t206;
    t508=t264.*t114;
    t509=t507+t508;
    t520=t247.*y1;
    t523=t279.*t154.*t206-t459.*t207-t439.*t427-0.2e1.*t264.*t520;
    t524=t275.*t523;
    
    e11=B1.*((t91-t99-t120-t121.*t163).*t166.*t168./0.8e1+t298.*t166.*t168./0.4e1)+B2.*(t338.*t166.*t168./0.8e1+t469.*t166.*t168./0.4e1)+B3.*(t475.*t163.*t166.*t168./0.8e1+(t173.*(-t484-t486+t489+t491)+t497-t226.*t504-t506.*t332.*t509.*y1-t506.*t335.*t509.*t160+t506.*t524).*t166.*t168./0.4e1);
    
    
    t533=0.1e1./y1.*t9;
    t536=0.1e1./t16.*t22;
    t545=(t131.*t33+t6.*t28.*t34-0.2e1.*t6.*t27.*t40.*t12).*t51;
    t547=0.1e1./t56.*t62;
    t549=t6.*t68;
    t556=(t151.*t72+t549.*t73-0.2e1.*t6.*t67.*t79.*t12).*t88;
    t558=t2.*(-0.2e1.*t533+t536+t545+t547+t556);
    t559=t103.*y2;
    t560=t108.*y2;
    t561=t112.*y2;
    t562=t117.*y2;
    t563=-t559-t560-t561-t562;
    t564=t100.*t563;
    t567=t136.*t129+t156.*t149;
    t569=t105.*t14;
    t570=y2.*t129;
    t574=t138.*y2;
    t576=t114.*t14;
    t577=y2.*t149;
    t581=t158.*y2;
    t583=t569.*t570-t133.*t570-t132.*t105.*t574+t576.*t577-t153.*t577-t152.*t114.*t581;
    t589=-t533+t547+t556;
    t590=t589.*t177;
    t593=t173.*t6;
    t594=t116.*t188;
    t597=y2.*t176;
    t598=t192.*t597;
    t599=y1.*t116;
    t600=t197.*y2;
    t602=t598+t599.*t600+t486;
    t605=t176.*t149;
    t606=t605.*t206;
    t607=t256.*t606;
    t608=t593.*t12;
    t609=t206.*t68;
    t614=a.*t6;
    t616=0.3e1.*t614.*t397;
    t617=t217.*t68;
    t619=t6.*t217;
    t624=t243.*y2;
    t625=t599.*t624;
    t626=t247.*y2;
    t627=t232.*t626;
    t635=t114.*t158;
    t644=t111.*y2;
    t645=t644.*t176;
    t647=t2.*t68;
    t648=t121.*t14;
    t658=t174.*t590+t362.*t188-t593.*t594.*t68+t180.*t96.*t602+t607-t608.*t205.*t609-t608.*t605.*t192+t395-t616+t617.*t235-t619.*t112.*t234-t619.*t117.*t234+t226.*t97.*(-t625-t486-0.2e1.*t627)+t617.*t268-t619.*t332.*t267-t619.*t635.*t267+t226.*t275.*(-t276.*t262.*t68.*y2+t253.*(t279.*y2.*t257.*t176+t283.*t645+t647.*t648)+0.2e1.*t290.*t291.*y2);
    t665=t69.*t96;
    t666=t29.*t129;
    t667=t69.*t149;
    t668=t666+t667;
    t670=t29.*t93+t665-t12.*t668;
    t673=t16.*t105;
    t676=t321.*y2;
    t679=t105.*t138.*y2;
    t681=t56.*t114;
    t684=t332.*y2;
    t686=t635.*y2;
    t698=t232.*t645;
    t700=t192.*y2;
    t710=t388.*y2;
    t718=t3.*t111;
    t724=t3.*t247;
    t736=t232.*t54;
    t737=t221.*t378;
    t741=y2.*t443;
    t745=t279.*t577.*t441;
    t748=t459.*t441.*t68.*y2;
    t749=t462.*t710;
    t762=t166.*t168;
    t769=t96.*t481;
    t770=t6.*t116;
    t772=t6.*t96;
    t773=t772.*t192;
    t774=t253.*t206;
    t775=t487.*t68;
    t776=t31.*t775;
    t777=t212.*t111;
    t778=t31.*t777;
    t781=t617.*t494;
    t783=t498.*y2;
    t785=t501.*y2;
    t787=t68.*(-0.2e1.*t783-t785);
    t789=t217.*t12;
    t790=t275.*t509;
    t792=t619.*t12;
    t799=t609.*y2;
    t804=t279.*t577.*t206-t459.*t799-t439.*t700-0.2e1.*t264.*t626;
    t805=t275.*t804;
    t814=-t322-t325-t333-t336;
    t822=t342-nu;
    t823=t822.*t68;
    t826=(t342+0.1e1-t172).*t12;
    t832=t193.*t183+t354-a+t772.*t186;
    t833=t832.*t68;
    t842=t173.*t56.*t176;
    t845=t402-a;
    t849=-t172+t68.*t845+t549.*t242+t614.*t111;
    t850=t849.*t68;
    t853=t111.*t845;
    t855=t6.*t111;
    t860=t6.*t247;
    t868=t217.*t158;
    t869=a.*t12;
    t870=t842+t869;
    t876=a.*t56;
    t877=t355.*t255;
    t879=t6.*t385-t876.*t877;
    t881=t385-t68.*t870+t264.*t56.*t176.*t111-t275.*t879;
    t884=t217.*t149;
    t885=t111.*t870;
    t890=t264.*t56;
    t898=t876.*t176;
    t899=t111.*t255;
    t900=t899.*y1;
    t902=t114.*t12;
    t903=t902.*y1;
    t909=t173.*(t823.*t185-t826.*t307)+t351.*t833.*y1-t362.*(t184+t366-t770.*t197.*y1-t772.*t427)-t607+t842.*t488+t842.*t214-t395+t399-t400.*t850.*y1+t412.*(-t853.*y1+t416-t855.*t231-t6.*t114.*t422.*y1-t860.*t425.*y1-0.3e1.*t614.*t221.*y1)-t868.*t881.*t160+t884.*(t885.*y1-t415.*t210+t264.*t210.*t111-0.3e1.*t890.*t223+t332.*t879.*y1+t335.*t879.*t160-t275.*(-t869.*t877+t898.*t900-t898.*t903));
    t915=t12.*t28;
    t917=t16.*t102;
    t919=t16.*t28;
    t922=t56.*t111;
    t924=t56.*t68;
    t932=t2.*t173;
    t936=-t183.*t176+t187;
    t943=t180.*t176;
    t944=t158.*t441;
    t947=t180.*t605;
    t949=t192.*t378.*y1;
    t951=nu.*y1;
    t953=0.2e1.*t951.*t96;
    t954=t68+t96;
    t957=t375-t953-t232.*t68.*t954;
    t958=t96.*t957;
    t961=t116.*t957;
    t965=0.2e1.*nu.*t96;
    t968=0.2e1.*nu.*t3.*t501;
    t970=t111.*t954;
    t978=t226.*t176;
    t982=t171.*t12+t442+t264.*t114.*t378;
    t989=t247.*t378;
    t996=t932.*t178-t181.*t936.*t68.*y1-t180.*t96.*t201+t943.*t944.*t160+t947.*t949+t225-t227.*t958.*y1-t238.*t961.*y1+t226.*t97.*(-t965+t968-t182.*t954+t246.*t970-t232.*t68.*(-t213-t502))-t978.*t332.*t982.*y1-t978.*t335.*t982.*t160+t978.*t275.*(t458-t461-t463-0.2e1.*t264.*t989.*y1);
    t1002=t173.*t14;
    t1004=t6.*t14;
    t1013=t56.*t158;
    t1015=t56.*t149;
    t1021=y1.*t217;
    t1025=0.1e1+t508;
    t1028=t68.*t255;
    t1030=t31.*t14-t876.*t1028;
    t1032=t14.*(t12-t182)+t924.*t1025-t275.*t1030;
    t1035=t14.*a;
    t1040=t56.*t221;
    
    e12=B1.*((t558-t312-t564-t12.*t567-t121.*t583).*t166.*t168./0.8e1+t658.*t166.*t168./0.4e1)./0.2e1+B2.*((t228.*t670+t3.*t563+t673.*t475.*t129-t320.*t676-t320.*t679+t681.*t475.*t149-t331.*t684-t331.*t686).*t166.*t168./0.8e1+(t173.*(t344.*t485-t347.*t667)-t351.*t359.*y2+t362.*(-t698-t196.*t600-t199.*t700)+t375.*t158.*t381.*t68.*y2-t375.*t149.*(-t245.*t475.*t378+t377.*t710)+t225-t400.*t409.*y2+t412.*(-t413.*y2+t718.*t242.*y2+t3.*t114.*t422.*y2+t724.*t425.*y2+0.3e1.*t246.*t221.*y2)-t396.*t158.*t445.*t68.*y2+t396.*t149.*(-0.3e1.*t736.*t737.*y2+t576.*t741-t153.*t741+t156.*(-t745+t748+t749))).*t166.*t168./0.4e1)./0.2e1+B3.*(t14.*t567.*t762./0.8e1+t475.*t583.*t166.*t168./0.8e1+(t173.*(t769-t770.*t482-t773-t774+t776+t778)-t781+t619.*t495-t226.*t787+t789.*t790-t792.*t332.*t509-t792.*t635.*t509+t506.*t805).*t166.*t168./0.4e1)./0.2e1+B1.*((t173.*t310-t6.*(-t104-t109-t113-t118-t12.*t814)).*t166.*t168./0.8e1+t909.*t166.*t168./0.4e1)./0.2e1+B2.*((t91+t99+t120-y2.*(t915.*t129-t917.*t134-t919.*t141+t279.*t149-t922.*t154-t924.*t161)).*t166.*t168./0.8e1+t996.*t166.*t168./0.4e1)./0.2e1+B3.*((t1002.*t308-t1004.*t814).*t166.*t168./0.8e1+(t173.*(-t14.*t160.*t149-t769+t196.*t482+t200+t774-t1013.*t207-t1015.*t427)+t781-t3.*t217.*t495+t1021.*t504+t868.*t1032.*t160-t884.*(t1035.*t213+t279.*t1025-t922.*t1025.*y1-0.2e1.*t1040.*t736+t332.*t1030.*y1+t335.*t1030.*t160-t275.*(-t869.*t1028+t876.*t900-t876.*t903))).*t166.*t168./0.4e1)./0.2e1;
    
    
    t1062=t19.*t14.*t22;
    t1068=(t29.*t34.*y3+t37.*t48.*y1).*t51;
    t1070=t59.*t14.*t62;
    t1073=0.2e1.*y3+(4.*a);
    t1080=(t69.*t73.*t1073./0.2e1-t76.*t85.*y1).*t88;
    t1082=t2.*(t1062+t1068-t1070+t1080);
    t1083=t103.*y3;
    t1085=t28.*y3;
    t1086=t1085-0.1e1;
    t1087=t28.*t107.*t1086;
    t1089=t112.*t1073./0.2e1;
    t1091=t68.*t1073./0.2e1;
    t1092=t1091+0.1e1;
    t1093=t501.*t1092;
    t1094=-t1083-t1087-t1089-t1093;
    t1095=t100.*t1094;
    t1096=y3.*t129;
    t1099=t1085-t12;
    t1100=t138.*t1099;
    t1102=t1073.*t149;
    t1107=t1091+t12;
    t1108=t158.*t1107;
    t1110=t569.*t1096-t133.*t1096-t136.*t1100+t576.*t1102./0.2e1-t153.*t1102./0.2e1-t156.*t1108;
    t1116=-t1070+t1080;
    t1117=t1116.*t177;
    t1124=t186.*t1092;
    t1126=t192.*t1073;
    t1128=t185.*t1126./0.2e1;
    t1129=t192.*t1073.*t176./0.2e1+t599.*t1124+t1128;
    t1132=t206.*t1107;
    t1135=t111.*t1073;
    t1136=t212.*t1135;
    t1139=t222.*t1073;
    t1141=0.3e1./0.2e1.*t218.*t1139;
    t1148=t230.*t1092;
    t1150=t247.*t1073;
    t1151=t232.*t1150;
    t1165=t279.*t1073./0.2e1+0.1e1;
    t1168=t1135.*t176;
    t1171=t14.*t1073;
    t1184=t174.*t1117-t180.*t594.*t1092+t180.*t96.*t1129-t204.*t205.*t1132-t211.*t1136./0.2e1+t598-t1141+t69.*t235-t227.*t235.*t1073./0.2e1-t271.*t239.*t1092+t226.*t97.*(-t599.*t1148-t1128-t1151)+t69.*t268-t227.*t268.*t1073./0.2e1-t271.*t272.*t1107+t226.*t275.*(-t276.*t262.*t1107+t253.*(t1165.*t257.*t176+t283.*t1168./0.2e1+t647.*t1171.*t12./0.2e1)-t869.*t176.*t114+t290.*t291.*t1073);
    t1191=t1092.*t96;
    t1192=t1099.*t129;
    t1193=t1107.*t149;
    t1194=t1192+t1193;
    t1196=t1086.*t93+t1191-t12.*t1194;
    t1203=t321.*y3;
    t1205=t324.*t1099;
    t1212=t332.*t1073;
    t1215=t335.*t1107;
    t1229=t232.*t1168./0.2e1;
    t1240=t388.*t1073;
    t1248=0.3e1./0.2e1.*t232.*t217.*t1139;
    t1254=t242.*t1073;
    t1257=t422.*t1092;
    t1259=t425.*t1073;
    t1262=t221.*t1073;
    t1273=t1073.*t443;
    t1278=t1165.*t149;
    t1279=t1278.*t441;
    t1281=t459.*t441.*t1107;
    t1283=t462.*t1240./0.2e1;
    t1289=t173.*(t343.*t1092.*t96-t347.*t1193)-t351.*t358.*t1092+t362.*(nu-t1229-t196.*t1124-t199.*t1126./0.2e1)+t375.*t382.*t1107-t375.*t149.*(t434-t245.*t1171.*t378./0.2e1+t377.*t1240./0.2e1)-t194+t1248+t96.*t408-t400.*t408.*t1092+t412.*(-t413.*t1073./0.2e1+t718.*t1254./0.2e1+t405.*t1257+t724.*t1259./0.2e1+0.3e1./0.2e1.*t246.*t1262)+t605.*t445-t396.*t446.*t1107+t396.*t149.*(t949-0.3e1./0.2e1.*t736.*t737.*t1073+t576.*t1273./0.2e1-t153.*t1273./0.2e1+t156.*(-t1279+t1281+t1283));
    t1299=t481.*t1092;
    t1302=t485.*t1126./0.2e1;
    t1303=t487.*t1107;
    t1304=t121.*t1303;
    t1306=t490.*t1126./0.2e1;
    t1310=t495.*t1073;
    t1313=t498.*t1073;
    t1314=t116.*t1092;
    t1316=t68.*(-t1313-t1314);
    t1330=t264.*t1150;
    t1331=t1278.*t206-t459.*t1132-t439.*t1126./0.2e1+t245-t1330;
    t1332=t275.*t1331;
    t1340=t102.*y1;
    t1341=t105.*t12;
    t1344=t27.*t12-y3;
    t1345=t1344.*t102;
    t1347=t1344.*t28;
    t1352=t1341.*t134-t1345.*t134-t1347.*t141-t902.*t154+t899.*t154+t1028.*t161;
    t1362=t965+t245;
    t1363=t111.*t1362;
    t1366=nu.*t116;
    t1372=0.1e1-t172-t507-t508;
    t1388=t228.*t14;
    t1397=t16.*t1344;
    t1401=t12.*t255.*t275;
    t1403=t681.*t13.*t149;
    t1404=t56.*t255;
    t1405=t1404.*t333;
    t1406=t1404.*t336;
    t1407=t1388.*(t306-t307)-t28+t68-y1.*(-t1340+t213)+t12.*t1344.*t316+t673.*t13.*t129-t1397.*t322-t1397.*t325-t1401-t1403+t1405+t1406;
    t1418=t2.*t3;
    t1423=t149.*t206;
    t1425=t2.*t56;
    t1427=t1425.*t149;
    t1429=t217.*t111;
    t1430=t375-t953-t233;
    t1435=t255.*t176;
    t1436=t438-t439;
    t1437=t68.*t1436;
    t1439=t54.*t56;
    t1442=t14-t1439.*t114-t1404.*t275;
    t1444=t434+t1435.*t1437+t182.*t1442;
    t1449=t111.*t1436;
    t1466=t174.*t176.*(t305-t12.*t160.*t149)-t2.*t96.*t230+t1418.*t422.*t68+t1418.*t425.*t111+t438.*t1423-t1425.*t488-t1427.*t427-t1429.*t1430.*y1+t617.*(-t965+t968-t245+t249)+t868.*t1444.*t160-t884.*(t902.*t193.*t1436-t1435.*t1449.*y1+t1435.*t68.*(-t279.*t154+t459.*t160)-t192.*t1442.*y1+t182.*(-t147.*t114+0.2e1.*t1439.*t520-t1401-t1403+t1405+t1406));
    t1481=t2.*y2.*t14;
    t1484=t226.*t14;
    t1485=0.1e1+t507+t508;
    
    e13=B1.*((t1082-t1095-t121.*t1110).*t166.*t168./0.8e1+t1184.*t166.*t168./0.4e1)./0.2e1+B2.*((t228.*t1196+t3.*t1094-t14.*t132.*t316+t673.*t15.*t129-t320.*t1203-t320.*t1205+t14.*t152.*t275+t681.*t1171.*t149./0.2e1-t331.*t1212./0.2e1-t331.*t1215).*t166.*t168./0.8e1+t1289.*t166.*t168./0.4e1)./0.2e1+B3.*(t475.*t1110.*t166.*t168./0.8e1+(t173.*(-t480.*t1299-t1302+t1304+t1306)-t69.*t494+t226.*t1310./0.2e1-t226.*t1316+t121.*t790-t506.*t332.*t509.*t1073./0.2e1-t506.*t335.*t509.*t1107+t506.*t1332).*t166.*t168./0.4e1)./0.2e1+B1.*(y2.*(-t1340+t213-t12.*t1352).*t762./0.8e1+(t2.*(t173.*t175.*t176-t625-t486+t489+t491)-t226.*t1363.*y1+t226.*t68.*(-0.2e1.*t1366.*t159-0.2e1.*t499)-t506.*t332.*t1372.*y1-t506.*t335.*t1372.*t160-t506.*t275.*t523).*t166.*t168./0.4e1)./0.2e1+B2.*(t1407.*t166.*t168./0.8e1+t1466.*t166.*t168./0.4e1)./0.2e1+B3.*((t2.*(-t24+t52+t64-t89)+t475.*t1352).*t166.*t168./0.8e1+(t2.*t175-t1481.*t488-t1481.*t214-t1484.*t332.*t1485.*y1-t1484.*t335.*t1485.*t160+t1484.*t524).*t166.*t168./0.4e1)./0.2e1;
    
    
    t1501=t316+t275;
    t1506=-t676-t679-t684-t686;
    t1522=t6.*y2;
    t1557=t879.*y2;
    t1562=t899.*y2;
    t1564=t902.*y2;
    t1590=t116.*t936;
    t1598=t593.*t176;
    t1608=t951.*t785;
    t1618=t275.*t982;
    t1620=t619.*t176;
    t1631=t932.*t590+t362.*t936-t593.*t1590.*t68-t180.*t96.*t602-t375.*t149.*t441+t1598.*t944.*t68+t1598.*t212.*t388-t395+t616+t617.*t958-t619.*t112.*t957-t619.*t117.*t957+t226.*t97.*(0.2e1.*t1608+t232.*t970.*y2-t232.*t68.*(-t644-t785))+t396.*t1618-t1620.*t332.*t982-t1620.*t635.*t982+t978.*t275.*(t745-t748-t749-0.2e1.*t264.*t989.*y2);
    t1660=t1030.*y2;
    
    e22=B1.*((t173.*t670-0.2e1.*y2.*(t94+t97-t12.*t1501)-t6.*(-t559-t560-t561-t562-t12.*t1506)).*t166.*t168./0.8e1+(t173.*(t823.*t485-t826.*t667)+t351.*t833.*y2-t362.*(t698+0.2e1.*t485.*t186-t1522.*t116.*t197-t1522.*t96.*t192)+t842.*t487.*t69+t842.*t212.*t644+t225-t400.*t850.*y2+t412.*(-t853.*y2+0.2e1.*t69.*t242-t1522.*t111.*t242-t1522.*t114.*t422-t1522.*t247.*t425+0.2e1.*t700-0.3e1.*a.*t1522.*t221)-t868.*t881.*t68.*y2+t884.*(t885.*y2-0.3e1.*t890.*t222.*y2+t332.*t1557+t635.*t1557-t275.*(0.2e1.*y2.*t385+t898.*t1562-t898.*t1564))).*t166.*t168./0.4e1)+B2.*((t558+t312+t564-t919.*t129-t924.*t149-y2.*(-t917.*t570-t673.*t574-t922.*t577-t681.*t581)).*t166.*t168./0.8e1+t1631.*t166.*t168./0.4e1)+B3.*((t1002.*t668-0.2e1.*t475.*t1501-t1004.*t1506).*t166.*t168./0.8e1+(t173.*(-t143.*t577+t484+t486-t1013.*t799-t1015.*t700)-t497+t1021.*t787+t868.*t1032.*t68.*y2-t884.*(t1035.*t644-t922.*t1025.*y2-0.2e1.*t1040.*t264.*y2+t332.*t1660+t635.*t1660-t275.*(0.2e1.*t648+t876.*t1562-t876.*t1564))).*t166.*t168./0.4e1);
    
    
    t1678=-t1203-t1205-t1212./0.2e1-t1215;
    t1721=t14.*t176;
    t1734=t899.*t1073;
    t1743=t173.*(t822.*t1092.*t96-t826.*t1193)+t351.*t832.*t1092-t362.*(t1229+nu-t770.*t1124-t772.*t1126./0.2e1)-t1002.*t606+t842.*t1303+t842.*t1136./0.2e1-t194+t1248+t96.*t849-t400.*t849.*t1092+t412.*(-t853.*t1073./0.2e1-t855.*t1254./0.2e1-t549.*t1257-t860.*t1259./0.2e1-0.3e1./0.2e1.*t614.*t1262)+t149.*t881-t868.*t881.*t1107+t884.*(t885.*t1073./0.2e1-t415.*t1721+t876.*t365+t264.*t1721.*t111-0.3e1./0.2e1.*t890.*t1139+t332.*t879.*t1073./0.2e1+t335.*t879.*t1107-t275.*(-t1035.*t877+t898.*t1734./0.2e1-t876.*t355.*t1165));
    t1781=0.2e1.*t951.*t1314;
    t1785=t1135./0.2e1;
    t1806=t932.*t1117-t180.*t1590.*t1092-t180.*t96.*t1129+t943.*t944.*t1107+t947.*t192.*t378.*t1073./0.2e1-t598+t1141+t69.*t958-t227.*t958.*t1073./0.2e1-t271.*t961.*t1092+t226.*t97.*(t1781+t232.*t970.*t1073./0.2e1-t232.*t68.*(-t1785-t1314))+t597.*t1618-t978.*t332.*t982.*t1073./0.2e1-t978.*t335.*t982.*t1107+t978.*t275.*(t1279-t1281-t1283+t245.*t378-t264.*t989.*t1073);
    t1851=t68.*t1165;
    t1865=t1347.*t129-t1028.*t149;
    t1871=t102.*y2;
    t1880=t1341.*t570-t1345.*t570-t1344.*t105.*t574-t902.*t577+t899.*t577+t255.*t114.*t581;
    t1898=t275.*t1372;
    t1921=t681.*t490;
    t1922=t1404.*t684;
    t1923=t1404.*t686;
    t1928=t279.*t577;
    t1932=t2.*y1;
    t1935=t1932.*t96;
    t1980=t2.*t14.*t1423;
    t1982=t2.*t6.*t14;
    t1986=t275.*t1485;
    t1988=t619.*t14;
    
    e23=B1.*((t173.*t1196-t6.*(-t1083-t1087-t1089-t1093-t12.*t1678)).*t166.*t168./0.8e1+t1743.*t166.*t168./0.4e1)./0.2e1+B2.*((t1082+t1095-y2.*(-t122.*t129-t917.*t1096-t919.*t1100+t143.*t149-t922.*t1102./0.2e1-t924.*t1108)).*t166.*t168./0.8e1+t1806.*t166.*t168./0.4e1)./0.2e1+B3.*((t1002.*t1194-t1004.*t1678).*t166.*t168./0.8e1+(t173.*(-t14.*t1107.*t149+t599.*t1299+t1128+t14.*t149.*t206-t1013.*t1132-t1015.*t1126./0.2e1)+t159.*t494-t1021.*t1310./0.2e1+t1021.*t1316-t149.*t1032+t868.*t1032.*t1107-t884.*(t1035.*t1135./0.2e1+t143.*t1025-t922.*t1025.*t1073./0.2e1+t924.*(t245-t1330)+t332.*t1030.*t1073./0.2e1+t335.*t1030.*t1107-t275.*(-t1035.*t1028+t876.*t1734./0.2e1-t876.*t1851))).*t166.*t168./0.4e1)./0.2e1+B1.*((t28-t68-t12.*t1865).*t166.*t168./0.8e1+y2.*(-t1871+t644-t12.*t1880).*t762./0.8e1+(t2.*(t173.*t589.*t176+t242-t770.*t243-t773-t774+t776+t778)+t617.*t1362-t619.*t1363+t226.*t68.*(-0.2e1.*t1366.*t69-0.2e1.*t783)+t789.*t1898-t792.*t332.*t1372-t792.*t635.*t1372-t506.*t275.*t804).*t166.*t168./0.4e1)./0.2e1+B2.*((t1388.*(t666-t667)-y1.*(-t1871+t644)+t673.*t121.*t129-t1397.*t676-t1397.*t679-t1921+t1922+t1923).*t166.*t168./0.8e1+(t174.*t176.*(t665-t1928)+t1932.*t116.*t624+t1935.*t700-t1425.*t158.*t799-t1427.*t700-t1429.*t1430.*y2+t617.*(0.2e1.*t1608+0.2e1.*t627)+t868.*t1444.*t68.*y2-t884.*(t902.*t597.*t1436-t1435.*t1449.*y2+t1435.*t68.*(-t1928+t459.*t69)-t192.*t1442.*y2+t182.*(0.2e1.*t1439.*t626-t1921+t1922+t1923))).*t166.*t168./0.4e1)./0.2e1+B3.*((t2.*(t536+t545-t547-t556)+t14.*t1865+t475.*t1880).*t166.*t168./0.8e1+(t2.*t589+t1980-t1982.*t775-t1982.*t777+t217.*t14.*t1986-t1988.*t332.*t1485-t1988.*t635.*t1485+t1484.*t805).*t166.*t168./0.4e1)./0.2e1;
    
    
    t2001=t102.*y3;
    t2003=t915.*y3-0.1e1;
    t2012=t2003.*t28.*t129-t1345.*t1096-t1347.*t1100-t1851.*t149+t899.*t1102./0.2e1+t1028.*t1108;
    t2060=t14.*t255.*t275;
    t2062=t56.*t1165.*t275;
    t2064=t1404.*t1212./0.2e1;
    t2065=t1404.*t1215;
    t2108=t174.*t176.*(t1191-t12.*t1107.*t149)+t1932.*t1257+t1935.*t1126./0.2e1+t1980-t1425.*t1303-t1427.*t1126./0.2e1+t68.*t1430-t1429.*t1430.*t1073./0.2e1+t617.*(t1781+t1151)-t149.*t1444+t868.*t1444.*t1107-t884.*(t1165.*t176.*t1437-t1435.*t1449.*t1073./0.2e1+t1435.*t68.*(-t1278+t459.*t1107)-t192.*t1442.*t1073./0.2e1+t182.*(-t681-t55.*t114+t1439.*t1150-t2060-t2062+t2064+t2065));
    
    e33=B1.*(y2.*(-t2001+t1785-t12.*t2012).*t762./0.8e1+(t2.*(t173.*t1116.*t176-t480.*t1148-t1302+t1304+t1306)+t69.*t1362-t226.*t1363.*t1073./0.2e1+t226.*t68.*(-0.2e1.*t1366.*t1092-t1313)+t121.*t1898-t506.*t332.*t1372.*t1073./0.2e1-t506.*t335.*t1372.*t1107-t506.*t275.*t1331).*t166.*t168./0.4e1)+B2.*((t1388.*(t1192-t1193)-y1.*(-t2001+t1785)-t14.*t1344.*t316+t16.*t2003.*t316-t1397.*t1203-t1397.*t1205-t2060-t2062+t2064+t2065).*t166.*t168./0.8e1+t2108.*t166.*t168./0.4e1)+B3.*((t2.*(t1062+t1068+t1070-t1080)+t475.*t2012).*t166.*t168./0.8e1+(t2.*t1116-t1481.*t1303-t1481.*t1136./0.2e1+t475.*t1986-t1484.*t332.*t1485.*t1073./0.2e1-t1484.*t335.*t1485.*t1107+t1484.*t1332).*t166.*t168./0.4e1);

end

