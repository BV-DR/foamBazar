#! /bin/bash
BC="constant/polyMesh/boundary"

if ! [ -f $BC ]; then
    echo $BC "doesn't exist."
    exit 1
fi

mkdir -p mk0.sh

# Create U
OUT="mk0.sh/U"; echo $OUT
    sed -e "s/polyBoundaryMesh/volVectorField/" \
        -e "s/constant\/polyMesh/0/" \
        -e "s/boundary/U/" \
        -e "s/(/{/" \
        -e "s/)/}/" $BC | \
    sed "17 i dimensions      [0 1 -1 0 0 0 0];
         18 i internalField   uniform (0 0 0);\n
         18 c boundaryField
         /nFaces/ c \        value           uniform (0 0 0);
         /startFace/ d" > $OUT
    sed -i -e "/top/,/}/s/patch;/pressureInletOutletVelocity;/" $OUT
    sed -i -e "s/wall;/fixedValue;/" $OUT

# Create pd
OUT="mk0.sh/pd"; echo $OUT
    sed -e "s/polyBoundaryMesh/volScalarField/" \
        -e "s/constant\/polyMesh/0/" \
        -e "s/boundary/pd/" \
        -e "s/(/{/" \
        -e "s/)/}/" $BC | \
    sed "17 i dimensions      [1 -1 -2 0 0 0 0];
         18 i internalField   uniform 0;\n
         18 c boundaryField
         /nFaces/ d;
         /startFace/ d" > $OUT
    sed -i -e "/top\|atmosphere\|Top\|TOP\|Atmosphere/,/}/s/patch;/totalPressure;\n\
        U               U;\n\
        phi             phi;\n\
        rho             none;\n\
        psi             none;\n\
        gamma           1;\n\
        p0              uniform 0;\n\
        value           uniform 0;/" $OUT
    sed -i -e "s/wall;/zeroGradient;/" $OUT

# Create p
OUT="mk0.sh/p"; echo $OUT
    cp "mk0.sh/pd" $OUT
    sed -i -e "s/object      pd;/object      p;/" $OUT
    sed -i -e "s/zeroGradient;/buoyantPressure;\n\
        gradient        uniform 0;\n\
        rho             rho;/" $OUT

# Create gamma
OUT="mk0.sh/gamma"; echo $OUT
    sed -e "s/polyBoundaryMesh/volScalarField/" \
        -e "s/constant\/polyMesh/0/" \
        -e "s/boundary/gamma/" \
        -e "s/(/{/" \
        -e "s/)/}/" $BC | \
    sed "17 i dimensions      [0 0 0 0 0 0 0];
         18 i internalField   uniform 0;\n
         18 c boundaryField
         /nFaces/ d;
         /startFace/ d" > $OUT
    sed -i -e "/top\|atmosphere\|Top\|TOP\|Atmosphere/,/}/s/patch;/inletOutlet;\n\
        inletValue      uniform 0;\n\
        value           uniform 0;/" $OUT
    sed -i -e "s/wall;/zeroGradient;/" $OUT

# Create alpha1
OUT="mk0.sh/alpha1"; echo $OUT
    cp "mk0.sh/gamma" $OUT
    sed -i -e "s/object      gamma;/object      alpha1;/" $OUT
    
# Create cellMotionUz
OUT="mk0.sh/cellMotionUz"; echo $OUT
    sed -e "s/polyBoundaryMesh/volScalarField/" \
        -e "s/constant\/polyMesh/0/" \
        -e "s/boundary/motionU/" \
        -e "s/(/{/" \
        -e "s/)/}/" $BC | \
    sed "17 i dimensions      [0 1 -1 0 0 0 0];
         18 i internalField   uniform 0;\n
         18 c boundaryField
         /nFaces/ d;
         /startFace/ d" > $OUT
    sed -i -e "/bottom\|top\|atmosphere\|Top\|TOP\|Atmosphere/,/}/s/patch\|wall;/fixedValue;\n\
        value           uniform 0;/" $OUT
    sed -i -e "/back\|front\|left\|right/,/}/s/wall;/slip;/" $OUT

# Create cellMotionU
OUT="mk0.sh/cellMotionU"; echo $OUT
    sed -e "s/polyBoundaryMesh/volVectorField/" \
        -e "s/constant\/polyMesh/0/" \
        -e "s/boundary/motionU/" \
        -e "s/(/{/" \
        -e "s/)/}/" $BC | \
    sed "17 i dimensions      [0 1 -1 0 0 0 0];
         18 i internalField   uniform (0 0 0);\n
         18 c boundaryField
         /nFaces/ d;
         /startFace/ d" > $OUT
    sed -i -e "/bottom\|top\|atmosphere\|Top\|TOP\|Atmosphere/,/}/s/patch\|wall;/fixedValue;\n\
        value           uniform (0 0 0);/" $OUT
    sed -i -e "/back\|front\|left\|right/,/}/s/wall;/slip;/" $OUT
    
# Create pointMotionUz
OUT="mk0.sh/pointMotionUz"; echo $OUT
	cp "mk0.sh/cellMotionUz" $OUT
    sed -i -e "s/volScalarField;/pointScalarField;/" $OUT
    sed -i -e "s/motionU;/pointMotionUz;/" $OUT 

# Create pointMotionU
OUT="mk0.sh/pointMotionU"; echo $OUT
    cp "mk0.sh/cellMotionU" $OUT
    sed -i -e "s/volVectorField;/pointVectorField;/" $OUT
    sed -i -e "s/motionU;/pointMotionU;/" $OUT           



