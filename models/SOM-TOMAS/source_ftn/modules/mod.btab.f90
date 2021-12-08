!     ==========================================================================
!       THIS IS THE MODULE FOR TABULATED COEFF. FOR DIFF.-LIMITED PARTITIONING
!     ==========================================================================

      MODULE hBTAB
  
      DOUBLE PRECISION L_TAB(61)
      DOUBLE PRECISION BETAS_TAB(61,10)

      DATA L_TAB(1) /-9.95e-01/ 
      DATA L_TAB(2) /-9.90e-01/ 
      DATA L_TAB(3) /-9.80e-01/ 
      DATA L_TAB(4) /-9.70e-01/ 
      DATA L_TAB(5) /-9.60e-01/ 
      DATA L_TAB(6) /-9.50e-01/ 
      DATA L_TAB(7) /-9.40e-01/ 
      DATA L_TAB(8) /-9.30e-01/ 
      DATA L_TAB(9) /-9.20e-01/ 
      DATA L_TAB(10) /-9.10e-01/ 
      DATA L_TAB(11) /-9.00e-01/ 
      DATA L_TAB(12) /-8.50e-01/ 
      DATA L_TAB(13) /-8.00e-01/ 
      DATA L_TAB(14) /-7.00e-01/ 
      DATA L_TAB(15) /-6.00e-01/ 
      DATA L_TAB(16) /-5.00e-01/ 
      DATA L_TAB(17) /-4.00e-01/ 
      DATA L_TAB(18) /-3.00e-01/ 
      DATA L_TAB(19) /-2.00e-01/ 
      DATA L_TAB(20) /-1.00e-01/ 
      DATA L_TAB(21) /0.00e+00/ 
      DATA L_TAB(22) /1.00e-01/ 
      DATA L_TAB(23) /2.00e-01/ 
      DATA L_TAB(24) /3.00e-01/ 
      DATA L_TAB(25) /4.00e-01/ 
      DATA L_TAB(26) /5.00e-01/ 
      DATA L_TAB(27) /6.00e-01/ 
      DATA L_TAB(28) /7.00e-01/ 
      DATA L_TAB(29) /8.00e-01/ 
      DATA L_TAB(30) /9.00e-01/ 
      DATA L_TAB(31) /1.00e+00/ 
      DATA L_TAB(32) /1.50e+00/ 
      DATA L_TAB(33) /2.00e+00/ 
      DATA L_TAB(34) /3.00e+00/ 
      DATA L_TAB(35) /4.00e+00/ 
      DATA L_TAB(36) /5.00e+00/ 
      DATA L_TAB(37) /6.00e+00/ 
      DATA L_TAB(38) /7.00e+00/ 
      DATA L_TAB(39) /8.00e+00/ 
      DATA L_TAB(40) /9.00e+00/ 
      DATA L_TAB(41) /1.00e+01/ 
      DATA L_TAB(42) /1.50e+01/ 
      DATA L_TAB(43) /2.00e+01/ 
      DATA L_TAB(44) /3.00e+01/ 
      DATA L_TAB(45) /4.00e+01/ 
      DATA L_TAB(46) /5.00e+01/ 
      DATA L_TAB(47) /6.00e+01/ 
      DATA L_TAB(48) /7.00e+01/ 
      DATA L_TAB(49) /8.00e+01/ 
      DATA L_TAB(50) /9.00e+01/ 
      DATA L_TAB(51) /1.00e+02/ 
      DATA L_TAB(52) /1.50e+02/ 
      DATA L_TAB(53) /2.00e+02/ 
      DATA L_TAB(54) /3.00e+02/ 
      DATA L_TAB(55) /4.00e+02/ 
      DATA L_TAB(56) /5.00e+02/ 
      DATA L_TAB(57) /6.00e+02/ 
      DATA L_TAB(58) /7.00e+02/ 
      DATA L_TAB(59) /8.00e+02/ 
      DATA L_TAB(60) /9.00e+02/ 
      DATA L_TAB(61) /1.00e+03/ 

      DATA BETAS_TAB(1,1) /0.12245/ 
      DATA BETAS_TAB(1,2) /4.49455/ 
      DATA BETAS_TAB(1,3) /7.72585/ 
      DATA BETAS_TAB(1,4) /10.90455/ 
      DATA BETAS_TAB(1,5) /14.06655/ 
      DATA BETAS_TAB(1,6) /17.22105/ 
      DATA BETAS_TAB(1,7) /20.37155/ 
      DATA BETAS_TAB(1,8) /23.51965/ 
      DATA BETAS_TAB(1,9) /26.66625/ 
      DATA BETAS_TAB(1,10) /29.81175/ 
      DATA BETAS_TAB(2,1) /0.17305/ 
      DATA BETAS_TAB(2,2) /4.49565/ 
      DATA BETAS_TAB(2,3) /7.72655/ 
      DATA BETAS_TAB(2,4) /10.90505/ 
      DATA BETAS_TAB(2,5) /14.06695/ 
      DATA BETAS_TAB(2,6) /17.22135/ 
      DATA BETAS_TAB(2,7) /20.37175/ 
      DATA BETAS_TAB(2,8) /23.51985/ 
      DATA BETAS_TAB(2,9) /26.66645/ 
      DATA BETAS_TAB(2,10) /29.81195/ 
      DATA BETAS_TAB(3,1) /0.24445/ 
      DATA BETAS_TAB(3,2) /4.49785/ 
      DATA BETAS_TAB(3,3) /7.72785/ 
      DATA BETAS_TAB(3,4) /10.90595/ 
      DATA BETAS_TAB(3,5) /14.06765/ 
      DATA BETAS_TAB(3,6) /17.22195/ 
      DATA BETAS_TAB(3,7) /20.37225/ 
      DATA BETAS_TAB(3,8) /23.52035/ 
      DATA BETAS_TAB(3,9) /26.66685/ 
      DATA BETAS_TAB(3,10) /29.81225/ 
      DATA BETAS_TAB(4,1) /0.29915/ 
      DATA BETAS_TAB(4,2) /4.50005/ 
      DATA BETAS_TAB(4,3) /7.72915/ 
      DATA BETAS_TAB(4,4) /10.90685/ 
      DATA BETAS_TAB(4,5) /14.06835/ 
      DATA BETAS_TAB(4,6) /17.22245/ 
      DATA BETAS_TAB(4,7) /20.37275/ 
      DATA BETAS_TAB(4,8) /23.52075/ 
      DATA BETAS_TAB(4,9) /26.66715/ 
      DATA BETAS_TAB(4,10) /29.81265/ 
      DATA BETAS_TAB(5,1) /0.34505/ 
      DATA BETAS_TAB(5,2) /4.50235/ 
      DATA BETAS_TAB(5,3) /7.73045/ 
      DATA BETAS_TAB(5,4) /10.90775/ 
      DATA BETAS_TAB(5,5) /14.06905/ 
      DATA BETAS_TAB(5,6) /17.22305/ 
      DATA BETAS_TAB(5,7) /20.37325/ 
      DATA BETAS_TAB(5,8) /23.52115/ 
      DATA BETAS_TAB(5,9) /26.66755/ 
      DATA BETAS_TAB(5,10) /29.81295/ 
      DATA BETAS_TAB(6,1) /0.38535/ 
      DATA BETAS_TAB(6,2) /4.50455/ 
      DATA BETAS_TAB(6,3) /7.73175/ 
      DATA BETAS_TAB(6,4) /10.90875/ 
      DATA BETAS_TAB(6,5) /14.06975/ 
      DATA BETAS_TAB(6,6) /17.22365/ 
      DATA BETAS_TAB(6,7) /20.37375/ 
      DATA BETAS_TAB(6,8) /23.52155/ 
      DATA BETAS_TAB(6,9) /26.66795/ 
      DATA BETAS_TAB(6,10) /29.81325/ 
      DATA BETAS_TAB(7,1) /0.42175/ 
      DATA BETAS_TAB(7,2) /4.50675/ 
      DATA BETAS_TAB(7,3) /7.73305/ 
      DATA BETAS_TAB(7,4) /10.90965/ 
      DATA BETAS_TAB(7,5) /14.07045/ 
      DATA BETAS_TAB(7,6) /17.22425/ 
      DATA BETAS_TAB(7,7) /20.37425/ 
      DATA BETAS_TAB(7,8) /23.52205/ 
      DATA BETAS_TAB(7,9) /26.66835/ 
      DATA BETAS_TAB(7,10) /29.81365/ 
      DATA BETAS_TAB(8,1) /0.45505/ 
      DATA BETAS_TAB(8,2) /4.50895/ 
      DATA BETAS_TAB(8,3) /7.73435/ 
      DATA BETAS_TAB(8,4) /10.91055/ 
      DATA BETAS_TAB(8,5) /14.07115/ 
      DATA BETAS_TAB(8,6) /17.22485/ 
      DATA BETAS_TAB(8,7) /20.37475/ 
      DATA BETAS_TAB(8,8) /23.52245/ 
      DATA BETAS_TAB(8,9) /26.66865/ 
      DATA BETAS_TAB(8,10) /29.81395/ 
      DATA BETAS_TAB(9,1) /0.48595/ 
      DATA BETAS_TAB(9,2) /4.51125/ 
      DATA BETAS_TAB(9,3) /7.73565/ 
      DATA BETAS_TAB(9,4) /10.91145/ 
      DATA BETAS_TAB(9,5) /14.07185/ 
      DATA BETAS_TAB(9,6) /17.22545/ 
      DATA BETAS_TAB(9,7) /20.37525/ 
      DATA BETAS_TAB(9,8) /23.52285/ 
      DATA BETAS_TAB(9,9) /26.66905/ 
      DATA BETAS_TAB(9,10) /29.81425/ 
      DATA BETAS_TAB(10,1) /0.51495/ 
      DATA BETAS_TAB(10,2) /4.51345/ 
      DATA BETAS_TAB(10,3) /7.73695/ 
      DATA BETAS_TAB(10,4) /10.91235/ 
      DATA BETAS_TAB(10,5) /14.07255/ 
      DATA BETAS_TAB(10,6) /17.22595/ 
      DATA BETAS_TAB(10,7) /20.37575/ 
      DATA BETAS_TAB(10,8) /23.52325/ 
      DATA BETAS_TAB(10,9) /26.66945/ 
      DATA BETAS_TAB(10,10) /29.81465/ 
      DATA BETAS_TAB(11,1) /0.54225/ 
      DATA BETAS_TAB(11,2) /4.51565/ 
      DATA BETAS_TAB(11,3) /7.73815/ 
      DATA BETAS_TAB(11,4) /10.91325/ 
      DATA BETAS_TAB(11,5) /14.07335/ 
      DATA BETAS_TAB(11,6) /17.22655/ 
      DATA BETAS_TAB(11,7) /20.37625/ 
      DATA BETAS_TAB(11,8) /23.52375/ 
      DATA BETAS_TAB(11,9) /26.66985/ 
      DATA BETAS_TAB(11,10) /29.81495/ 
      DATA BETAS_TAB(12,1) /0.66085/ 
      DATA BETAS_TAB(12,2) /4.52675/ 
      DATA BETAS_TAB(12,3) /7.74465/ 
      DATA BETAS_TAB(12,4) /10.91785/ 
      DATA BETAS_TAB(12,5) /14.07685/ 
      DATA BETAS_TAB(12,6) /17.22945/ 
      DATA BETAS_TAB(12,7) /20.37865/ 
      DATA BETAS_TAB(12,8) /23.52585/ 
      DATA BETAS_TAB(12,9) /26.67165/ 
      DATA BETAS_TAB(12,10) /29.81665/ 
      DATA BETAS_TAB(13,1) /0.75935/ 
      DATA BETAS_TAB(13,2) /4.53785/ 
      DATA BETAS_TAB(13,3) /7.75115/ 
      DATA BETAS_TAB(13,4) /10.92245/ 
      DATA BETAS_TAB(13,5) /14.08045/ 
      DATA BETAS_TAB(13,6) /17.23235/ 
      DATA BETAS_TAB(13,7) /20.38115/ 
      DATA BETAS_TAB(13,8) /23.52795/ 
      DATA BETAS_TAB(13,9) /26.67355/ 
      DATA BETAS_TAB(13,10) /29.81835/ 
      DATA BETAS_TAB(14,1) /0.92075/ 
      DATA BETAS_TAB(14,2) /4.56005/ 
      DATA BETAS_TAB(14,3) /7.76405/ 
      DATA BETAS_TAB(14,4) /10.93165/ 
      DATA BETAS_TAB(14,5) /14.08755/ 
      DATA BETAS_TAB(14,6) /17.23815/ 
      DATA BETAS_TAB(14,7) /20.38605/ 
      DATA BETAS_TAB(14,8) /23.53225/ 
      DATA BETAS_TAB(14,9) /26.67735/ 
      DATA BETAS_TAB(14,10) /29.82165/ 
      DATA BETAS_TAB(15,1) /1.05275/ 
      DATA BETAS_TAB(15,2) /4.58215/ 
      DATA BETAS_TAB(15,3) /7.77695/ 
      DATA BETAS_TAB(15,4) /10.94075/ 
      DATA BETAS_TAB(15,5) /14.09465/ 
      DATA BETAS_TAB(15,6) /17.24395/ 
      DATA BETAS_TAB(15,7) /20.39095/ 
      DATA BETAS_TAB(15,8) /23.53645/ 
      DATA BETAS_TAB(15,9) /26.68105/ 
      DATA BETAS_TAB(15,10) /29.82505/ 
      DATA BETAS_TAB(16,1) /1.16555/ 
      DATA BETAS_TAB(16,2) /4.60425/ 
      DATA BETAS_TAB(16,3) /7.78985/ 
      DATA BETAS_TAB(16,4) /10.94995/ 
      DATA BETAS_TAB(16,5) /14.10175/ 
      DATA BETAS_TAB(16,6) /17.24975/ 
      DATA BETAS_TAB(16,7) /20.39585/ 
      DATA BETAS_TAB(16,8) /23.54075/ 
      DATA BETAS_TAB(16,9) /26.68485/ 
      DATA BETAS_TAB(16,10) /29.82835/ 
      DATA BETAS_TAB(17,1) /1.26445/ 
      DATA BETAS_TAB(17,2) /4.62615/ 
      DATA BETAS_TAB(17,3) /7.80275/ 
      DATA BETAS_TAB(17,4) /10.95905/ 
      DATA BETAS_TAB(17,5) /14.10885/ 
      DATA BETAS_TAB(17,6) /17.25555/ 
      DATA BETAS_TAB(17,7) /20.40075/ 
      DATA BETAS_TAB(17,8) /23.54495/ 
      DATA BETAS_TAB(17,9) /26.68855/ 
      DATA BETAS_TAB(17,10) /29.83175/ 
      DATA BETAS_TAB(18,1) /1.35255/ 
      DATA BETAS_TAB(18,2) /4.64795/ 
      DATA BETAS_TAB(18,3) /7.81565/ 
      DATA BETAS_TAB(18,4) /10.96825/ 
      DATA BETAS_TAB(18,5) /14.11595/ 
      DATA BETAS_TAB(18,6) /17.26135/ 
      DATA BETAS_TAB(18,7) /20.40565/ 
      DATA BETAS_TAB(18,8) /23.54925/ 
      DATA BETAS_TAB(18,9) /26.69225/ 
      DATA BETAS_TAB(18,10) /29.83505/ 
      DATA BETAS_TAB(19,1) /1.43205/ 
      DATA BETAS_TAB(19,2) /4.66955/ 
      DATA BETAS_TAB(19,3) /7.82845/ 
      DATA BETAS_TAB(19,4) /10.97735/ 
      DATA BETAS_TAB(19,5) /14.12305/ 
      DATA BETAS_TAB(19,6) /17.26715/ 
      DATA BETAS_TAB(19,7) /20.41055/ 
      DATA BETAS_TAB(19,8) /23.55345/ 
      DATA BETAS_TAB(19,9) /26.69605/ 
      DATA BETAS_TAB(19,10) /29.83845/ 
      DATA BETAS_TAB(20,1) /1.50445/ 
      DATA BETAS_TAB(20,2) /4.69105/ 
      DATA BETAS_TAB(20,3) /7.84125/ 
      DATA BETAS_TAB(20,4) /10.98645/ 
      DATA BETAS_TAB(20,5) /14.13005/ 
      DATA BETAS_TAB(20,6) /17.27295/ 
      DATA BETAS_TAB(20,7) /20.41545/ 
      DATA BETAS_TAB(20,8) /23.55775/ 
      DATA BETAS_TAB(20,9) /26.69975/ 
      DATA BETAS_TAB(20,10) /29.84175/ 
      DATA BETAS_TAB(21,1) /1.57075/ 
      DATA BETAS_TAB(21,2) /4.71235/ 
      DATA BETAS_TAB(21,3) /7.85395/ 
      DATA BETAS_TAB(21,4) /10.99555/ 
      DATA BETAS_TAB(21,5) /14.13715/ 
      DATA BETAS_TAB(21,6) /17.27875/ 
      DATA BETAS_TAB(21,7) /20.42035/ 
      DATA BETAS_TAB(21,8) /23.56195/ 
      DATA BETAS_TAB(21,9) /26.70355/ 
      DATA BETAS_TAB(21,10) /29.84515/ 
      DATA BETAS_TAB(22,1) /1.63195/ 
      DATA BETAS_TAB(22,2) /4.73355/ 
      DATA BETAS_TAB(22,3) /7.86665/ 
      DATA BETAS_TAB(22,4) /11.00465/ 
      DATA BETAS_TAB(22,5) /14.14425/ 
      DATA BETAS_TAB(22,6) /17.28455/ 
      DATA BETAS_TAB(22,7) /20.42525/ 
      DATA BETAS_TAB(22,8) /23.56615/ 
      DATA BETAS_TAB(22,9) /26.70725/ 
      DATA BETAS_TAB(22,10) /29.84845/ 
      DATA BETAS_TAB(23,1) /1.68865/ 
      DATA BETAS_TAB(23,2) /4.75445/ 
      DATA BETAS_TAB(23,3) /7.87935/ 
      DATA BETAS_TAB(23,4) /11.01375/ 
      DATA BETAS_TAB(23,5) /14.15125/ 
      DATA BETAS_TAB(23,6) /17.29035/ 
      DATA BETAS_TAB(23,7) /20.43015/ 
      DATA BETAS_TAB(23,8) /23.57045/ 
      DATA BETAS_TAB(23,9) /26.71105/ 
      DATA BETAS_TAB(23,10) /29.85185/ 
      DATA BETAS_TAB(24,1) /1.74135/ 
      DATA BETAS_TAB(24,2) /4.77515/ 
      DATA BETAS_TAB(24,3) /7.89195/ 
      DATA BETAS_TAB(24,4) /11.02275/ 
      DATA BETAS_TAB(24,5) /14.15835/ 
      DATA BETAS_TAB(24,6) /17.29615/ 
      DATA BETAS_TAB(24,7) /20.43505/ 
      DATA BETAS_TAB(24,8) /23.57465/ 
      DATA BETAS_TAB(24,9) /26.71475/ 
      DATA BETAS_TAB(24,10) /29.85515/ 
      DATA BETAS_TAB(25,1) /1.79055/ 
      DATA BETAS_TAB(25,2) /4.79565/ 
      DATA BETAS_TAB(25,3) /7.90455/ 
      DATA BETAS_TAB(25,4) /11.03185/ 
      DATA BETAS_TAB(25,5) /14.16535/ 
      DATA BETAS_TAB(25,6) /17.30185/ 
      DATA BETAS_TAB(25,7) /20.43995/ 
      DATA BETAS_TAB(25,8) /23.57895/ 
      DATA BETAS_TAB(25,9) /26.71855/ 
      DATA BETAS_TAB(25,10) /29.85855/ 
      DATA BETAS_TAB(26,1) /1.83655/ 
      DATA BETAS_TAB(26,2) /4.81585/ 
      DATA BETAS_TAB(26,3) /7.91705/ 
      DATA BETAS_TAB(26,4) /11.04085/ 
      DATA BETAS_TAB(26,5) /14.17245/ 
      DATA BETAS_TAB(26,6) /17.30765/ 
      DATA BETAS_TAB(26,7) /20.44485/ 
      DATA BETAS_TAB(26,8) /23.58315/ 
      DATA BETAS_TAB(26,9) /26.72225/ 
      DATA BETAS_TAB(26,10) /29.86185/ 
      DATA BETAS_TAB(27,1) /1.87975/ 
      DATA BETAS_TAB(27,2) /4.83585/ 
      DATA BETAS_TAB(27,3) /7.92955/ 
      DATA BETAS_TAB(27,4) /11.04985/ 
      DATA BETAS_TAB(27,5) /14.17945/ 
      DATA BETAS_TAB(27,6) /17.31345/ 
      DATA BETAS_TAB(27,7) /20.44965/ 
      DATA BETAS_TAB(27,8) /23.58735/ 
      DATA BETAS_TAB(27,9) /26.72595/ 
      DATA BETAS_TAB(27,10) /29.86525/ 
      DATA BETAS_TAB(28,1) /1.92035/ 
      DATA BETAS_TAB(28,2) /4.85555/ 
      DATA BETAS_TAB(28,3) /7.94185/ 
      DATA BETAS_TAB(28,4) /11.05875/ 
      DATA BETAS_TAB(28,5) /14.18645/ 
      DATA BETAS_TAB(28,6) /17.31915/ 
      DATA BETAS_TAB(28,7) /20.45455/ 
      DATA BETAS_TAB(28,8) /23.59165/ 
      DATA BETAS_TAB(28,9) /26.72975/ 
      DATA BETAS_TAB(28,10) /29.86855/ 
      DATA BETAS_TAB(29,1) /1.95855/ 
      DATA BETAS_TAB(29,2) /4.87505/ 
      DATA BETAS_TAB(29,3) /7.95425/ 
      DATA BETAS_TAB(29,4) /11.06775/ 
      DATA BETAS_TAB(29,5) /14.19345/ 
      DATA BETAS_TAB(29,6) /17.32495/ 
      DATA BETAS_TAB(29,7) /20.45945/ 
      DATA BETAS_TAB(29,8) /23.59585/ 
      DATA BETAS_TAB(29,9) /26.73345/ 
      DATA BETAS_TAB(29,10) /29.87195/ 
      DATA BETAS_TAB(30,1) /1.99465/ 
      DATA BETAS_TAB(30,2) /4.89425/ 
      DATA BETAS_TAB(30,3) /7.96645/ 
      DATA BETAS_TAB(30,4) /11.07665/ 
      DATA BETAS_TAB(30,5) /14.20045/ 
      DATA BETAS_TAB(30,6) /17.33065/ 
      DATA BETAS_TAB(30,7) /20.46435/ 
      DATA BETAS_TAB(30,8) /23.60005/ 
      DATA BETAS_TAB(30,9) /26.73715/ 
      DATA BETAS_TAB(30,10) /29.87525/ 
      DATA BETAS_TAB(31,1) /2.02875/ 
      DATA BETAS_TAB(31,2) /4.91315/ 
      DATA BETAS_TAB(31,3) /7.97865/ 
      DATA BETAS_TAB(31,4) /11.08555/ 
      DATA BETAS_TAB(31,5) /14.20745/ 
      DATA BETAS_TAB(31,6) /17.33635/ 
      DATA BETAS_TAB(31,7) /20.46915/ 
      DATA BETAS_TAB(31,8) /23.60425/ 
      DATA BETAS_TAB(31,9) /26.74095/ 
      DATA BETAS_TAB(31,10) /29.87855/ 
      DATA BETAS_TAB(32,1) /2.17465/ 
      DATA BETAS_TAB(32,2) /5.00365/ 
      DATA BETAS_TAB(32,3) /8.03845/ 
      DATA BETAS_TAB(32,4) /11.12955/ 
      DATA BETAS_TAB(32,5) /14.24215/ 
      DATA BETAS_TAB(32,6) /17.36495/ 
      DATA BETAS_TAB(32,7) /20.49345/ 
      DATA BETAS_TAB(32,8) /23.62535/ 
      DATA BETAS_TAB(32,9) /26.75955/ 
      DATA BETAS_TAB(32,10) /29.89525/ 
      DATA BETAS_TAB(33,1) /2.28895/ 
      DATA BETAS_TAB(33,2) /5.08695/ 
      DATA BETAS_TAB(33,3) /8.09615/ 
      DATA BETAS_TAB(33,4) /11.17275/ 
      DATA BETAS_TAB(33,5) /14.27635/ 
      DATA BETAS_TAB(33,6) /17.39325/ 
      DATA BETAS_TAB(33,7) /20.51755/ 
      DATA BETAS_TAB(33,8) /23.64635/ 
      DATA BETAS_TAB(33,9) /26.77805/ 
      DATA BETAS_TAB(33,10) /29.91185/ 
      DATA BETAS_TAB(34,1) /2.45565/ 
      DATA BETAS_TAB(34,2) /5.23295/ 
      DATA BETAS_TAB(34,3) /8.20455/ 
      DATA BETAS_TAB(34,4) /11.25605/ 
      DATA BETAS_TAB(34,5) /14.34335/ 
      DATA BETAS_TAB(34,6) /17.44905/ 
      DATA BETAS_TAB(34,7) /20.56525/ 
      DATA BETAS_TAB(34,8) /23.68795/ 
      DATA BETAS_TAB(34,9) /26.81495/ 
      DATA BETAS_TAB(34,10) /29.94495/ 
      DATA BETAS_TAB(35,1) /2.57045/ 
      DATA BETAS_TAB(35,2) /5.35405/ 
      DATA BETAS_TAB(35,3) /8.30295/ 
      DATA BETAS_TAB(35,4) /11.33485/ 
      DATA BETAS_TAB(35,5) /14.40795/ 
      DATA BETAS_TAB(35,6) /17.50345/ 
      DATA BETAS_TAB(35,7) /20.61205/ 
      DATA BETAS_TAB(35,8) /23.72895/ 
      DATA BETAS_TAB(35,9) /26.85145/ 
      DATA BETAS_TAB(35,10) /29.97775/ 
      DATA BETAS_TAB(36,1) /2.65365/ 
      DATA BETAS_TAB(36,2) /5.45435/ 
      DATA BETAS_TAB(36,3) /8.39135/ 
      DATA BETAS_TAB(36,4) /11.40865/ 
      DATA BETAS_TAB(36,5) /14.46985/ 
      DATA BETAS_TAB(36,6) /17.55625/ 
      DATA BETAS_TAB(36,7) /20.65785/ 
      DATA BETAS_TAB(36,8) /23.76925/ 
      DATA BETAS_TAB(36,9) /26.88735/ 
      DATA BETAS_TAB(36,10) /30.01025/ 
      DATA BETAS_TAB(37,1) /2.71645/ 
      DATA BETAS_TAB(37,2) /5.53785/ 
      DATA BETAS_TAB(37,3) /8.47025/ 
      DATA BETAS_TAB(37,4) /11.47725/ 
      DATA BETAS_TAB(37,5) /14.52885/ 
      DATA BETAS_TAB(37,6) /17.60715/ 
      DATA BETAS_TAB(37,7) /20.70245/ 
      DATA BETAS_TAB(37,8) /23.80885/ 
      DATA BETAS_TAB(37,9) /26.92285/ 
      DATA BETAS_TAB(37,10) /30.04225/ 
      DATA BETAS_TAB(38,1) /2.76535/ 
      DATA BETAS_TAB(38,2) /5.60775/ 
      DATA BETAS_TAB(38,3) /8.54055/ 
      DATA BETAS_TAB(38,4) /11.54075/ 
      DATA BETAS_TAB(38,5) /14.58465/ 
      DATA BETAS_TAB(38,6) /17.65625/ 
      DATA BETAS_TAB(38,7) /20.74575/ 
      DATA BETAS_TAB(38,8) /23.84745/ 
      DATA BETAS_TAB(38,9) /26.95755/ 
      DATA BETAS_TAB(38,10) /30.07385/ 
      DATA BETAS_TAB(39,1) /2.80445/ 
      DATA BETAS_TAB(39,2) /5.66685/ 
      DATA BETAS_TAB(39,3) /8.60305/ 
      DATA BETAS_TAB(39,4) /11.59935/ 
      DATA BETAS_TAB(39,5) /14.63735/ 
      DATA BETAS_TAB(39,6) /17.70315/ 
      DATA BETAS_TAB(39,7) /20.78775/ 
      DATA BETAS_TAB(39,8) /23.88515/ 
      DATA BETAS_TAB(39,9) /26.99165/ 
      DATA BETAS_TAB(39,10) /30.10485/ 
      DATA BETAS_TAB(40,1) /2.83635/ 
      DATA BETAS_TAB(40,2) /5.71725/ 
      DATA BETAS_TAB(40,3) /8.65875/ 
      DATA BETAS_TAB(40,4) /11.65325/ 
      DATA BETAS_TAB(40,5) /14.68695/ 
      DATA BETAS_TAB(40,6) /17.74805/ 
      DATA BETAS_TAB(40,7) /20.82825/ 
      DATA BETAS_TAB(40,8) /23.92175/ 
      DATA BETAS_TAB(40,9) /27.02505/ 
      DATA BETAS_TAB(40,10) /30.13535/ 
      DATA BETAS_TAB(41,1) /2.86275/ 
      DATA BETAS_TAB(41,2) /5.76055/ 
      DATA BETAS_TAB(41,3) /8.70835/ 
      DATA BETAS_TAB(41,4) /11.70265/ 
      DATA BETAS_TAB(41,5) /14.73345/ 
      DATA BETAS_TAB(41,6) /17.79085/ 
      DATA BETAS_TAB(41,7) /20.86725/ 
      DATA BETAS_TAB(41,8) /23.95735/ 
      DATA BETAS_TAB(41,9) /27.05755/ 
      DATA BETAS_TAB(41,10) /30.16525/ 
      DATA BETAS_TAB(42,1) /2.94755/ 
      DATA BETAS_TAB(42,2) /5.90795/ 
      DATA BETAS_TAB(42,3) /8.88975/ 
      DATA BETAS_TAB(42,4) /11.89585/ 
      DATA BETAS_TAB(42,5) /14.92505/ 
      DATA BETAS_TAB(42,6) /17.97425/ 
      DATA BETAS_TAB(42,7) /21.03975/ 
      DATA BETAS_TAB(42,8) /24.11835/ 
      DATA BETAS_TAB(42,9) /27.20735/ 
      DATA BETAS_TAB(42,10) /30.30475/ 
      DATA BETAS_TAB(43,1) /2.99305/ 
      DATA BETAS_TAB(43,2) /5.99205/ 
      DATA BETAS_TAB(43,3) /9.00185/ 
      DATA BETAS_TAB(43,4) /12.02505/ 
      DATA BETAS_TAB(43,5) /15.06245/ 
      DATA BETAS_TAB(43,6) /18.11365/ 
      DATA BETAS_TAB(43,7) /21.17715/ 
      DATA BETAS_TAB(43,8) /24.25155/ 
      DATA BETAS_TAB(43,9) /27.33515/ 
      DATA BETAS_TAB(43,10) /30.42665/ 
      DATA BETAS_TAB(44,1) /3.04055/ 
      DATA BETAS_TAB(44,2) /6.08315/ 
      DATA BETAS_TAB(44,3) /9.12935/ 
      DATA BETAS_TAB(44,4) /12.18065/ 
      DATA BETAS_TAB(44,5) /15.23795/ 
      DATA BETAS_TAB(44,6) /18.30175/ 
      DATA BETAS_TAB(44,7) /21.37215/ 
      DATA BETAS_TAB(44,8) /24.44895/ 
      DATA BETAS_TAB(44,9) /27.53185/ 
      DATA BETAS_TAB(44,10) /30.62025/ 
      DATA BETAS_TAB(45,1) /3.06515/ 
      DATA BETAS_TAB(45,2) /6.13105/ 
      DATA BETAS_TAB(45,3) /9.19875/ 
      DATA BETAS_TAB(45,4) /12.26875/ 
      DATA BETAS_TAB(45,5) /15.34175/ 
      DATA BETAS_TAB(45,6) /18.41805/ 
      DATA BETAS_TAB(45,7) /21.49795/ 
      DATA BETAS_TAB(45,8) /24.58165/ 
      DATA BETAS_TAB(45,9) /27.66915/ 
      DATA BETAS_TAB(45,10) /30.76035/ 
      DATA BETAS_TAB(46,1) /3.08005/ 
      DATA BETAS_TAB(46,2) /6.16055/ 
      DATA BETAS_TAB(46,3) /9.24205/ 
      DATA BETAS_TAB(46,4) /12.32465/ 
      DATA BETAS_TAB(46,5) /15.40905/ 
      DATA BETAS_TAB(46,6) /18.49525/ 
      DATA BETAS_TAB(46,7) /21.58365/ 
      DATA BETAS_TAB(46,8) /24.67435/ 
      DATA BETAS_TAB(46,9) /27.76735/ 
      DATA BETAS_TAB(46,10) /30.86295/ 
      DATA BETAS_TAB(47,1) /3.09015/ 
      DATA BETAS_TAB(47,2) /6.18055/ 
      DATA BETAS_TAB(47,3) /9.27145/ 
      DATA BETAS_TAB(47,4) /12.36315/ 
      DATA BETAS_TAB(47,5) /15.45585/ 
      DATA BETAS_TAB(47,6) /18.54975/ 
      DATA BETAS_TAB(47,7) /21.64495/ 
      DATA BETAS_TAB(47,8) /24.74165/ 
      DATA BETAS_TAB(47,9) /27.83995/ 
      DATA BETAS_TAB(47,10) /30.93985/ 
      DATA BETAS_TAB(48,1) /3.09735/ 
      DATA BETAS_TAB(48,2) /6.19495/ 
      DATA BETAS_TAB(48,3) /9.29275/ 
      DATA BETAS_TAB(48,4) /12.39115/ 
      DATA BETAS_TAB(48,5) /15.49015/ 
      DATA BETAS_TAB(48,6) /18.58995/ 
      DATA BETAS_TAB(48,7) /21.69065/ 
      DATA BETAS_TAB(48,8) /24.79235/ 
      DATA BETAS_TAB(48,9) /27.89515/ 
      DATA BETAS_TAB(48,10) /30.99905/ 
      DATA BETAS_TAB(49,1) /3.10285/ 
      DATA BETAS_TAB(49,2) /6.20575/ 
      DATA BETAS_TAB(49,3) /9.30895/ 
      DATA BETAS_TAB(49,4) /12.41245/ 
      DATA BETAS_TAB(49,5) /15.51635/ 
      DATA BETAS_TAB(49,6) /18.62085/ 
      DATA BETAS_TAB(49,7) /21.72595/ 
      DATA BETAS_TAB(49,8) /24.83175/ 
      DATA BETAS_TAB(49,9) /27.93835/ 
      DATA BETAS_TAB(49,10) /31.04575/ 
      DATA BETAS_TAB(50,1) /3.10705/ 
      DATA BETAS_TAB(50,2) /6.21425/ 
      DATA BETAS_TAB(50,3) /9.32155/ 
      DATA BETAS_TAB(50,4) /12.42915/ 
      DATA BETAS_TAB(50,5) /15.53705/ 
      DATA BETAS_TAB(50,6) /18.64525/ 
      DATA BETAS_TAB(50,7) /21.75395/ 
      DATA BETAS_TAB(50,8) /24.86325/ 
      DATA BETAS_TAB(50,9) /27.97295/ 
      DATA BETAS_TAB(50,10) /31.08335/ 
      DATA BETAS_TAB(51,1) /3.11045/ 
      DATA BETAS_TAB(51,2) /6.22105/ 
      DATA BETAS_TAB(51,3) /9.33175/ 
      DATA BETAS_TAB(51,4) /12.44255/ 
      DATA BETAS_TAB(51,5) /15.55365/ 
      DATA BETAS_TAB(51,6) /18.66505/ 
      DATA BETAS_TAB(51,7) /21.77675/ 
      DATA BETAS_TAB(51,8) /24.88885/ 
      DATA BETAS_TAB(51,9) /28.00135/ 
      DATA BETAS_TAB(51,10) /31.11425/ 
      DATA BETAS_TAB(52,1) /3.12075/ 
      DATA BETAS_TAB(52,2) /6.24155/ 
      DATA BETAS_TAB(52,3) /9.36245/ 
      DATA BETAS_TAB(52,4) /12.48335/ 
      DATA BETAS_TAB(52,5) /15.60435/ 
      DATA BETAS_TAB(52,6) /18.72535/ 
      DATA BETAS_TAB(52,7) /21.84655/ 
      DATA BETAS_TAB(52,8) /24.96785/ 
      DATA BETAS_TAB(52,9) /28.08925/ 
      DATA BETAS_TAB(52,10) /31.21075/ 
      DATA BETAS_TAB(53,1) /3.12595/ 
      DATA BETAS_TAB(53,2) /6.25195/ 
      DATA BETAS_TAB(53,3) /9.37795/ 
      DATA BETAS_TAB(53,4) /12.50395/ 
      DATA BETAS_TAB(53,5) /15.62995/ 
      DATA BETAS_TAB(53,6) /18.75605/ 
      DATA BETAS_TAB(53,7) /21.88215/ 
      DATA BETAS_TAB(53,8) /25.00835/ 
      DATA BETAS_TAB(53,9) /28.13455/ 
      DATA BETAS_TAB(53,10) /31.26085/ 
      DATA BETAS_TAB(54,1) /3.13115/ 
      DATA BETAS_TAB(54,2) /6.26235/ 
      DATA BETAS_TAB(54,3) /9.39345/ 
      DATA BETAS_TAB(54,4) /12.52465/ 
      DATA BETAS_TAB(54,5) /15.65585/ 
      DATA BETAS_TAB(54,6) /18.78705/ 
      DATA BETAS_TAB(54,7) /21.91825/ 
      DATA BETAS_TAB(54,8) /25.04945/ 
      DATA BETAS_TAB(54,9) /28.18065/ 
      DATA BETAS_TAB(54,10) /31.31195/ 
      DATA BETAS_TAB(55,1) /3.13375/ 
      DATA BETAS_TAB(55,2) /6.26755/ 
      DATA BETAS_TAB(55,3) /9.40125/ 
      DATA BETAS_TAB(55,4) /12.53505/ 
      DATA BETAS_TAB(55,5) /15.66885/ 
      DATA BETAS_TAB(55,6) /18.80255/ 
      DATA BETAS_TAB(55,7) /21.93635/ 
      DATA BETAS_TAB(55,8) /25.07015/ 
      DATA BETAS_TAB(55,9) /28.20395/ 
      DATA BETAS_TAB(55,10) /31.33775/ 
      DATA BETAS_TAB(56,1) /3.13535/ 
      DATA BETAS_TAB(56,2) /6.27065/ 
      DATA BETAS_TAB(56,3) /9.40595/ 
      DATA BETAS_TAB(56,4) /12.54125/ 
      DATA BETAS_TAB(56,5) /15.67665/ 
      DATA BETAS_TAB(56,6) /18.81195/ 
      DATA BETAS_TAB(56,7) /21.94725/ 
      DATA BETAS_TAB(56,8) /25.08265/ 
      DATA BETAS_TAB(56,9) /28.21795/ 
      DATA BETAS_TAB(56,10) /31.35335/ 
      DATA BETAS_TAB(57,1) /3.13635/ 
      DATA BETAS_TAB(57,2) /6.27275/ 
      DATA BETAS_TAB(57,3) /9.40905/ 
      DATA BETAS_TAB(57,4) /12.54545/ 
      DATA BETAS_TAB(57,5) /15.68185/ 
      DATA BETAS_TAB(57,6) /18.81825/ 
      DATA BETAS_TAB(57,7) /21.95455/ 
      DATA BETAS_TAB(57,8) /25.09095/ 
      DATA BETAS_TAB(57,9) /28.22735/ 
      DATA BETAS_TAB(57,10) /31.36375/ 
      DATA BETAS_TAB(58,1) /3.13715/ 
      DATA BETAS_TAB(58,2) /6.27425/ 
      DATA BETAS_TAB(58,3) /9.41135/ 
      DATA BETAS_TAB(58,4) /12.54845/ 
      DATA BETAS_TAB(58,5) /15.68555/ 
      DATA BETAS_TAB(58,6) /18.82265/ 
      DATA BETAS_TAB(58,7) /21.95975/ 
      DATA BETAS_TAB(58,8) /25.09695/ 
      DATA BETAS_TAB(58,9) /28.23405/ 
      DATA BETAS_TAB(58,10) /31.37115/ 
      DATA BETAS_TAB(59,1) /3.13765/ 
      DATA BETAS_TAB(59,2) /6.27535/ 
      DATA BETAS_TAB(59,3) /9.41305/ 
      DATA BETAS_TAB(59,4) /12.55065/ 
      DATA BETAS_TAB(59,5) /15.68835/ 
      DATA BETAS_TAB(59,6) /18.82605/ 
      DATA BETAS_TAB(59,7) /21.96375/ 
      DATA BETAS_TAB(59,8) /25.10135/ 
      DATA BETAS_TAB(59,9) /28.23905/ 
      DATA BETAS_TAB(59,10) /31.37675/ 
      DATA BETAS_TAB(60,1) /3.13815/ 
      DATA BETAS_TAB(60,2) /6.27625/ 
      DATA BETAS_TAB(60,3) /9.41435/ 
      DATA BETAS_TAB(60,4) /12.55245/ 
      DATA BETAS_TAB(60,5) /15.69055/ 
      DATA BETAS_TAB(60,6) /18.82865/ 
      DATA BETAS_TAB(60,7) /21.96675/ 
      DATA BETAS_TAB(60,8) /25.10485/ 
      DATA BETAS_TAB(60,9) /28.24295/ 
      DATA BETAS_TAB(60,10) /31.38105/ 
      DATA BETAS_TAB(61,1) /3.13845/ 
      DATA BETAS_TAB(61,2) /6.27695/ 
      DATA BETAS_TAB(61,3) /9.41535/ 
      DATA BETAS_TAB(61,4) /12.55385/ 
      DATA BETAS_TAB(61,5) /15.69225/ 
      DATA BETAS_TAB(61,6) /18.83075/ 
      DATA BETAS_TAB(61,7) /21.96915/ 
      DATA BETAS_TAB(61,8) /25.10765/ 
      DATA BETAS_TAB(61,9) /28.24605/ 
      DATA BETAS_TAB(61,10) /31.38455/ 

      END MODULE hBTAB