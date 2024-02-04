
particledia = sort(particledia, rev=true)
L1 = dimen1 #mm
L2 = dimen2 #mm
L3 = dimen3 #mm
Lmin = da/2+da/2+0.2*d0 #mm diameter of subcubes from division equals the possible largest distance between particle centers 
Nl1 = Int64(round(L1 / Lmin-0.5))
Nl2 = Int64(round(L2 / Lmin-0.5))
Nl3 = Int64(round(L3 / Lmin-0.5))
Li = maximum([L1 / Nl1, L2 / Nl2, L3 / Nl3]) 
magnifyp = LDPM.geometry_parameters[12] #control the miminum particle distance which equals manifyp*(radius of the first particle + radius of the second particle)
#thirdindexL = Int(round(Li/magnifyp/d0+0.5))
#thirdindexS = Int(round(Li/(magnifyp*d0)+0.5)^2*2)
thirdindexC = Int(round(Li/(d0/2+d0/2+0.2*d0)+0.5)^3*3)
point = [0, 0, 0]
pointsfinal = Vector{Float64}[]
pointsdiameterfinal= Float64[]
pointsinorder = Array{Any}(missing, Nl1, Nl2, Nl3, thirdindexC)
diametersinorder = Array{Any}(missing, Nl1, Nl2, Nl3, thirdindexC)
# eight corners
pointsinorder[1, 1, 1, 1] = [0.0, 0.0, 0.0] #
pointsinorder[Nl1, 1, 1, 1] = [L1, 0.0, 0.0]
pointsinorder[1, Nl2, 1, 1] = [0.0, L2, 0.0]
pointsinorder[Nl1, Nl2, 1, 1] = [L1, L2, 0.0]
pointsinorder[1, 1, Nl3, 1] = [0.0, 0.0, L3]
pointsinorder[Nl1, 1, Nl3, 1] = [L1, 0.0, L3]
pointsinorder[1, Nl2, Nl3, 1] = [0.0, L2, L3]
pointsinorder[Nl1, Nl2, Nl3, 1] = [L1, L2, L3]
diametersinorder[1, 1, 1, 1] = magnifyp*d0#
diametersinorder[Nl1, 1, 1, 1] = magnifyp*d0
diametersinorder[1, Nl2, 1, 1] = magnifyp*d0
diametersinorder[Nl1, Nl2, 1, 1] = magnifyp*d0
diametersinorder[1, 1, Nl3, 1] = magnifyp*d0
diametersinorder[Nl1, 1, Nl3, 1] = magnifyp*d0
diametersinorder[1, Nl2, Nl3, 1] = magnifyp*d0
diametersinorder[Nl1, Nl2, Nl3, 1] = magnifyp*d0
#random points on 12 edges
#x,0,0 edge
for i in 1:10000
    point = L1 * [rand(), 0, 0]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for ll in 1:thirdindexC
            if pointsinorder[ii, 1, 1, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[ii, 1, 1, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, 1, 1, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, 1, 1, mm] = point
                diametersinorder[boxnumx, 1, 1, mm] = magnifyp*d0
                break
            end
        end
    end
end
# #x,0,L edge
for i in 1:10000
    point =  [L1 *rand(), 0, L3]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for ll in 1:thirdindexC
            if pointsinorder[ii, 1, Nl3, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[ii, 1, Nl3, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, 1, Nl3, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, 1, Nl3, mm] = point
                diametersinorder[boxnumx, 1, Nl3, mm] = magnifyp*d0
                break
            end
        end
    end
    
end
# #x,L,0 edge
for i in 1:10000
    point = [L1 *rand(), L2, 0]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for ll in 1:thirdindexC
            if pointsinorder[ii, Nl2, 1, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[ii, Nl2, 1, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, Nl2, 1, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, Nl2, 1, mm] = point
                diametersinorder[boxnumx, Nl2, 1, mm] = magnifyp*d0
                break
            end
        end
    end
end
# #x,L,L edge
for i in 1:10000
    point =  [L1 *rand(), L2, L3]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for ll in 1:thirdindexC
            if pointsinorder[ii, Nl2, Nl3, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[ii, Nl2, Nl3, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, Nl2, Nl3, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, Nl2, Nl3, mm] = point
                diametersinorder[boxnumx, Nl2, Nl3, mm] = magnifyp*d0
                break
            end
        end
    end
end
#0,y,0 edge
for i in 1:10000
    point = L2 * [0, rand(), 0]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for ll in 1:thirdindexC
            if pointsinorder[1, jj, 1, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[1, jj, 1, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[1, boxnumy, 1, mm] isa Vector{Float64}

            else
                pointsinorder[1, boxnumy, 1, mm] = point
                diametersinorder[1, boxnumy, 1, mm] = magnifyp*d0
                break
            end
        end
    end
end
# #0,y,L edge
for i in 1:10000
    point =  [0, L2*rand(), L3]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for ll in 1:thirdindexC
            if pointsinorder[1, jj, Nl3, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[1, jj, Nl3, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[1, boxnumy, Nl3, mm] isa Vector{Float64}

            else
                pointsinorder[1, boxnumy, Nl3, mm] = point
                diametersinorder[1, boxnumy, Nl3, mm] = magnifyp*d0
                break
            end
        end
    end
    
end
# #L,y,0 edge
for i in 1:10000
    point = [L1, L2 * rand(), 0]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for ll in 1:thirdindexC
            if pointsinorder[Nl1, jj, 1, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[Nl1, jj, 1, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[Nl1, boxnumy, 1, mm] isa Vector{Float64}

            else
                pointsinorder[Nl1, boxnumy, 1, mm] = point
                diametersinorder[Nl1, boxnumy, 1, mm] = magnifyp*d0
                break
            end
        end
    end
    
end
# #L,y,L edge
for i in 1:10000
    point = [L1, L2 * rand(), L3]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for ll in 1:thirdindexC
            if pointsinorder[Nl1, jj, Nl3, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[Nl1, jj, Nl3, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[Nl1, boxnumy, Nl3, mm] isa Vector{Float64}

            else
                pointsinorder[Nl1, boxnumy, Nl3, mm] = point
                diametersinorder[Nl1, boxnumy, Nl3, mm] = magnifyp*d0
                break
            end
        end
    end
    
end
#0,0,z edge
for i in 1:10000
    point = L3 * [0, 0, rand()]
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end
    control = 1  #the points is effective
    for kk in kkinitial:kkend
        for ll in 1:thirdindexC
            if pointsinorder[1, 1, kk, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[1, 1, kk, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[1, 1, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[1, 1, boxnumz, mm] = point
                diametersinorder[1, 1, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end
end
#L,0,z edge
for i in 1:10000
    point = [L1, 0, L3 * rand()]
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end
    control = 1  #the points is effective
    for kk in kkinitial:kkend
        for ll in 1:thirdindexC
            if pointsinorder[Nl1, 1, kk, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[Nl1, 1, kk, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[Nl1, 1, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[Nl1, 1, boxnumz, mm] = point
                diametersinorder[Nl1, 1, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end
end
#0,L,z edge
for i in 1:10000
    point = [0, L2, L3 * rand()]
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end
    control = 1  #the points is effective
    for kk in kkinitial:kkend
        for ll in 1:thirdindexC
            if pointsinorder[1, Nl2, kk, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[1, Nl2, kk, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[1, Nl2, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[1, Nl2, boxnumz, mm] = point
                diametersinorder[1, Nl2, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end
end
#L,L,z edge
for i in 1:10000
    point = [L1, L2, L3 * rand()]
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end
    control = 1  #the points is effective
    for kk in kkinitial:kkend
        for ll in 1:thirdindexC
            if pointsinorder[Nl1, Nl2, kk, ll] isa Vector{Float64}
                Dis = euclidean(point, pointsinorder[Nl1, Nl2, kk, ll])
                if Dis < magnifyp*d0
                    control = 0 #the point should be deleted
                end

            end


        end

    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[Nl1, Nl2, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[Nl1, Nl2, boxnumz, mm] = point
                diametersinorder[Nl1, Nl2, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end
end
# random points in 6 surfaces 
#in surface xy
for i in 1:100000

    boxnumx = Int64(round(point[1] / Li + 0.5))
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end

    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for jj in jjinitial:jjend
            for ll in 1:thirdindexC
                if pointsinorder[ii, jj, 1, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[ii, jj, 1, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end


            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, boxnumy, 1, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, boxnumy, 1, mm] = point
                diametersinorder[boxnumx, boxnumy, 1, mm] = magnifyp*d0
                break
            end
        end
    end

end
for i in 1:100000

    boxnumx = Int64(round(point[1] / Li + 0.5))
    boxnumy = Int64(round(point[2] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end

    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for jj in jjinitial:jjend
            for ll in 1:thirdindexC
                if pointsinorder[ii, jj, Nl3, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[ii, jj, Nl3, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end


            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, boxnumy, Nl3, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, boxnumy, Nl3, mm] = point
                diametersinorder[boxnumx, boxnumy, Nl3, mm] = magnifyp*d0
                break
            end
        end
    end

end
# random points in surfaces yz
for i in 1:100000
    point = [0, L2 * rand(), L3 * rand()]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end

    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end

    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for kk in kkinitial:kkend
            for ll in 1:thirdindexC
                if pointsinorder[1, jj, kk, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[1, jj, kk, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end

            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[1, boxnumy, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[1, boxnumy, boxnumz, mm] = point
                diametersinorder[1, boxnumy, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end

end
for i in 1:100000
    point = [L1, L2 * rand(), L3 * rand()]
    boxnumy = Int64(round(point[2] / Li + 0.5))
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end

    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end

    control = 1  #the points is effective
    for jj in jjinitial:jjend
        for kk in kkinitial:kkend
            for ll in 1:thirdindexC
                if pointsinorder[Nl1, jj, kk, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[Nl1, jj, kk, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end


            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[Nl1, boxnumy, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[Nl1, boxnumy, boxnumz, mm] = point
                diametersinorder[Nl1, boxnumy, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end

end
# random points in surfaces xz
for i in 1:100000
    point = [L1 * rand(), 0, L3 * rand()]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end

    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end

    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for kk in kkinitial:kkend
            for ll in 1:thirdindexC
                if pointsinorder[ii, 1, kk, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[ii, 1, kk, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end


            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, 1, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, 1, boxnumz, mm] = point
                diametersinorder[boxnumx, 1, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end

end
for i in 1:100000
    point = [L1 * rand(), L2, L3 * rand()]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end

    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end

    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for kk in kkinitial:kkend
            for ll in 1:thirdindexC
                if pointsinorder[ii, Nl2, kk, ll] isa Vector{Float64}
                    Dis = euclidean(point, pointsinorder[ii, Nl2, kk, ll])
                    if Dis < magnifyp*d0
                        control = 0 #the point should be deleted
                    end

                end


            end

        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, Nl2, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, Nl2, boxnumz, mm] = point
                diametersinorder[boxnumx, Nl2, boxnumz, mm] = magnifyp*d0
                break
            end
        end
    end

end
#point in cubic
gdg = 0
while gdg<length(particledia)
    point = [L1 * rand(), L2 * rand(), L3 * rand()]
    boxnumx = Int64(round(point[1] / Li + 0.5))
    boxnumy = Int64(round(point[2] / Li + 0.5))
    boxnumz = Int64(round(point[3] / Li + 0.5))
    if boxnumx == 1
        iiinitial = 1
    else
        iiinitial = boxnumx - 1
    end
    if boxnumy == 1
        jjinitial = 1
    else
        jjinitial = boxnumy - 1
    end
    if boxnumz == 1
        kkinitial = 1
    else
        kkinitial = boxnumz - 1
    end
    if boxnumx == Nl1
        iiend = Nl1
    else
        iiend = boxnumx + 1
    end
    if boxnumy == Nl2
        jjend = Nl2
    else
        jjend = boxnumy + 1
    end
    if boxnumz == Nl3
        kkend = Nl3
    else
        kkend = boxnumz + 1
    end
    control = 1  #the points is effective
    for ii in iiinitial:iiend
        for jj in jjinitial:jjend
            for kk in kkinitial:kkend
                for ll in 1:thirdindexC
                    if pointsinorder[ii, jj, kk, ll] isa Vector{Float64}
                        Dis = euclidean(point, pointsinorder[ii, jj, kk, ll])
                        if Dis < (particledia[gdg+1]/2+diametersinorder[ii, jj, kk, ll]/2)*magnifyp
                            control = 0 #the point should be deleted
                        end

                    end


                end
            end
        end
    end
    if control == 1
        for mm in 1:thirdindexC
            if pointsinorder[boxnumx, boxnumy, boxnumz, mm] isa Vector{Float64}

            else
                pointsinorder[boxnumx, boxnumy, boxnumz, mm] = point
                diametersinorder[boxnumx, boxnumy, boxnumz, mm] = particledia[gdg+1]
                gdg = gdg+1
                break
            end
        end
    end
println(gdg)
end
#transfer positions of points
oo = 1
for iii in 1:Nl1
    for jjj in 1:Nl2
        for kkk in 1:Nl3
            for lll in 1:thirdindexC
                if pointsinorder[iii, jjj, kkk, lll] isa Vector{Float64}
                    push!(pointsfinal, pointsinorder[iii, jjj, kkk, lll])
                    oo = oo + 1
                end
            end
        end
    end
end
oo = 1
for iii in 1:Nl1
    for jjj in 1:Nl2
        for kkk in 1:Nl3
            for lll in 1:thirdindexC
                if diametersinorder[iii, jjj, kkk, lll] isa Float64
                    push!(pointsdiameterfinal, diametersinorder[iii, jjj, kkk, lll])
                    oo = oo + 1
                end
            end
        end
    end
end

# pointsfinal = pointsfinal[numbering1]
# pointsdiameterfinal = pointsdiameterfinal[numbering1]
pointsresult = mapreduce(permutedims, vcat, pointsfinal) # transfer Array to matrix
for klk = 1:length(pointsfinal)
    if pointsfinal[klk][1]==0 || pointsfinal[klk][2]==0 || pointsfinal[klk][3]==0||pointsfinal[klk][1]==L1||pointsfinal[klk][2]==L2||pointsfinal[klk][3]==L3
        pointsdiameterfinal[klk] = 0
    end
end

using JLD2
@save "geometry for concrete with steel_for show steel elements2.jld2" pointsresult pointsdiameterfinal
#Length(pointsfinal)
# using PlotlyJS
# function sphere1(r, C)   # r: radius; C: center [cx,cy,cz]
#     global n = 40
#     u = range(-π, π; length = n)
#     v = range(0, π; length = n)
#     x = C[1] .+ r*cos.(u) * sin.(v)'
#     y = C[2] .+ r*sin.(u) * sin.(v)'
#     z = C[3] .+ r*ones(n) * cos.(v)'
#     return x, y, z
# end
##pointsradiusfinal = pointsdiameterfinal / 2
# traceball = GenericTrace{Dict{Symbol,Any}}[]
# for defg = 1:Int(length(pointsfinal))
#     if pointsradiusfinal[defg] != 0
#     kok = sphere1(pointsradiusfinal[defg], pointsfinal[defg])#pointsradiusfinal[defg], pointsfinal[defg])
#     push!(traceball,PlotlyJS.surface(x=kok[1],y=kok[2],z=kok[3],showscale=false,surfacecolor=zeros(n,n)))
#     end
# end
# traceframe1 = PlotlyJS.scatter3d(; x=[0, 0], y=[0, 0], z=[0,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe2 = PlotlyJS.scatter3d(; x=[0, 0], y=[0, L], z=[0,0], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe3 = PlotlyJS.scatter3d(; x=[0, L], y=[0, 0], z=[0,0], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe4 = PlotlyJS.scatter3d(; x=[L, L], y=[L, L], z=[0,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe5 = PlotlyJS.scatter3d(; x=[L, L], y=[0, L], z=[L,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe6 = PlotlyJS.scatter3d(; x=[0, L], y=[L, L], z=[L,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe7 = PlotlyJS.scatter3d(; x=[0, 0], y=[0, L], z=[L,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe8 = PlotlyJS.scatter3d(; x=[0, L], y=[0, 0], z=[L,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe9 = PlotlyJS.scatter3d(; x=[L, L], y=[0, 0], z=[0,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe10 = PlotlyJS.scatter3d(; x=[L, L], y=[0, L], z=[0,0], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe11 = PlotlyJS.scatter3d(; x=[0, 0], y=[L, L], z=[0,L], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# traceframe12 = PlotlyJS.scatter3d(; x=[0, L], y=[L, L], z=[0,0], showlegend=false, mode="lines",marker=attr(color="rgb(127, 127, 127)"))
# layout = Layout(autosize=false, width=500, height=500,
#     margin=attr(l=40, r=20, b=20, t=65), legend=attr(
#         x=1.2,
#         y=1,
#         yanchor="bottom",
#         xanchor="right",
#         orientation="h"
#     ), scene_camera = attr(
#         up=attr(x=0, y=0, z=1),
#         center=attr(x=0, y=0, z=0),
#         eye=attr(x=1.7, y=1.2, z=1.2)))
#particle_position = plot([traceball; traceframe1; traceframe2;traceframe3;traceframe4;traceframe5;traceframe6;traceframe7;traceframe8;traceframe9;traceframe10;traceframe11;traceframe12], layout)
        #savefig(particle_position, "particle position with shape.svg")

# using PlotlyJS
#plt3d= PlotlyJS.plot3D(pointsresult[:,1],pointsresult[:,2], pointsresult[:,3], seriestype=:scatter, markersize = 7)
# plt3d = PlotlyJS.scatter3d(; x=pointsresult[:, 1], y=pointsresult[:, 2], z=pointsresult[:, 3],text=collect(1:length(pointsresult)/3),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))
# layout = Layout(margin=attr(l=40, r=40, t=40, b=40), scene_camera = attr(
#     up=attr(x=0, y=0, z=1),
#     center=attr(x=0, y=0, z=0),
#     eye=attr(x=1.7, y=1.2, z=1.2)))
# points_position = plot(plt3d, layout)
# # F̄[sequence_on_bottom_first]
# pointsresult4 = pointsresult[pointsresult[:, 3] .== 268.0, :]
# #pointsresult4 = pointsresult4[pointsresult4[:, 2] .!= 25.0, :]
# top_controll_points = pointsresult4[[9,20,31,36],:]
# top_controll_points = [top_controll_points[i,:] for i = 1:4]
# plt3d = PlotlyJS.scatter3d(; x=pointsresult4[:, 1], y=pointsresult4[:, 2], z=pointsresult4[:, 3],text=collect(1:length(pointsresult4)),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))
# layout = Layout(margin=attr(l=40, r=40, t=40, b=40), scene_camera = attr(
#     up=attr(x=0, y=0, z=1),
#     center=attr(x=0, y=0, z=0),
#     eye=attr(x=1.7, y=1.2, z=1.2)))
# points_position = plot(plt3d, layout)
# savefig(points_position, "points position.svg")
#rm("randompoints_mechanical_150.xls")
# using DelimitedFiles
# writedlm("randompoints_Kupher.xls", pointsresult)
#u----plot([plt3d; traceball; traceframe1; traceframe2;traceframe3;traceframe4;traceframe5;traceframe6;traceframe7;traceframe8;traceframe9;traceframe10;traceframe11;traceframe12], layout)#u-------------

#run(`paraview /path/to/your/file.vtk`)