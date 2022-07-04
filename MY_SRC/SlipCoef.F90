FUNCTION CalculSlc_haf (model, nodenumber, VarIn) RESULT(VarOut)
 USE DefUtils
 IMPLICIT NONE 
 TYPE(Model_t) :: model
 INTEGER :: nodenumber
 REAL (KIND=dp) :: VarIn(3) ! slipcoef, zs, bedrock
 REAL (KIND=dp) :: VarOut   ! slc

 TYPE(ValueList_t), POINTER :: material
 REAL(kind=dp) :: lda
 REAL(kind=dp) :: h_af
 REAL(kind=dp) :: h_af_init
 REAL(kind=dp) :: hth, zsl, rhoi, rhow
 
 ! inquire SSA friction exponent from Material properties
 material => GetMaterial()
 IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('SlipCoef_USF', "No material found?")
 ENDIF

 ! get needed variable
 hth = ListGetConstReal( Material, 'SSA Friction Threshold Height',UnFoundFatal=.TRUE.) 
 zsl = ListGetCReal( Model % Constants, 'Sea Level',UnFoundFatal=.TRUE.)
 rhoi = ListGetCReal( Model % Constants, 'Ice density',UnFoundFatal=.TRUE.)
 rhow = ListGetCReal( Model % Constants, 'water density',UnFoundFatal=.TRUE.)

 ! compute heigh above flotation (h_af)
 ! to distinguish later on flating and grounded part, haf use bedrock and not zb
 !    freeboard = ((rhow/rhoi)*(zsl - bedrock)+bedrock) 
 !    h_af = zs - freeboard
 h_af = VarIn(2) - ((rhow/rhoi)*(zsl - VarIn(3))+VarIn(3))

 ! PM : h_af_init commented as not used 
 !h_af = VarIn(2) - ((rhow/rhoi)*(zsl - VarIn(4))+VarIn(4))
 !h_af_init = MAX(0.1, VarIn(3) - ((rhow/rhoi)*(zsl - VarIn(4))+VarIn(4)))
 ! end PM

 !IF (h_af.GT.75.) THEN
 ! lda=+1.0
 !ELSE
 ! VarOut=h_af/MIN(75.0,h_af_init)
 !ENDIF

 ! compute scale scale factor  
 !     0 on floating part
 !     1 where above treshold
 ! [0-1] where between free floating and treshold
 !lda=h_af/MIN(h_af_init,hth)

 ! general case
 lda=h_af/hth

 ! where floating 
 IF (lda.LE.0.) THEN
 lda=+0.0
 ENDIF

 ! haf superior to h_threshold 
 IF (h_af.GE.hth) THEN
 lda=+1.0
 ENDIF

 ! output
 ! transition between both case)
 VarOut=VarIn(1)*lda
 
 END FUNCTION CalculSlc_haf
