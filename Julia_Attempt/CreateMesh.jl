using Gmsh

function build_mesh(escape_w, escape_b, fin_gap, n_fins_override=nothing;
                    outfile="coldplate.msh", nz=10,
                    dy_hi=3.0)   # ← coarse Y spacing in tall fluid region (mm)

    fin_t      = 0.2
    fin_pitch  = fin_t + fin_gap
    W_total    = 58.0
    H_fluid    = 32.8
    H_fin      = 3.0
    L_fin      = 53.0

    array_span = W_total - escape_w - escape_b
    n_fins     = isnothing(n_fins_override) ? floor(Int, array_span / fin_pitch) : n_fins_override

    println("n_fins=$(n_fins)  pitch=$(round(fin_pitch,digits=3))mm  array_span=$(round(array_span,digits=2))mm")

    x0 = 0.0;  x1 = escape_b;  x2 = W_total - escape_w;  x3 = W_total
    y0 = 0.0;  y1 = H_fin;     y2 = H_fluid
    z0 = 0.0;  z1 = L_fin

    println("Layout: escape_b=$(x1)mm | porous=$(x1) to $(x2)mm | escape_w=$(x3-x2)mm")

    gmsh.initialize()
    gmsh.model.add("coldplate")
    gmsh.option.setNumber("General.Terminal", 0)

    t_left_lo  = gmsh.model.occ.addBox(x0, y0, z0, x1-x0, y1-y0, z1)
    t_left_hi  = gmsh.model.occ.addBox(x0, y1, z0, x1-x0, y2-y1, z1)
    t_porous   = gmsh.model.occ.addBox(x1, y0, z0, x2-x1, y1-y0, z1)
    t_above    = gmsh.model.occ.addBox(x1, y1, z0, x2-x1, y2-y1, z1)
    t_right_lo = gmsh.model.occ.addBox(x2, y0, z0, x3-x2, y1-y0, z1)
    t_right_hi = gmsh.model.occ.addBox(x2, y1, z0, x3-x2, y2-y1, z1)

    all_tags = [t_left_lo, t_left_hi, t_porous, t_above, t_right_lo, t_right_hi]
    gmsh.model.occ.fragment([(3,t) for t in all_tags], [])
    gmsh.model.occ.synchronize()

    function find_vol(cx, cy, cz; tol=0.5)
        for (dim, tag) in gmsh.model.getEntities(3)
            xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
            if abs((xmin+xmax)/2-cx)<tol && abs((ymin+ymax)/2-cy)<tol && abs((zmin+zmax)/2-cz)<tol
                return tag
            end
        end
        error("Volume not found near centroid ($cx, $cy, $cz)")
    end

    cz = z1/2
    b_left_lo  = find_vol((x0+x1)/2, (y0+y1)/2, cz)
    b_left_hi  = find_vol((x0+x1)/2, (y1+y2)/2, cz)
    b_porous   = find_vol((x1+x2)/2, (y0+y1)/2, cz)
    b_above    = find_vol((x1+x2)/2, (y1+y2)/2, cz)
    b_right_lo = find_vol((x2+x3)/2, (y0+y1)/2, cz)
    b_right_hi = find_vol((x2+x3)/2, (y1+y2)/2, cz)

    all_boxes = [b_left_lo, b_left_hi, b_porous, b_above, b_right_lo, b_right_hi]
    names     = ["left_lo","left_hi","porous","above","right_lo","right_hi"]

    print("Volume tags —")
    for (n,t) in zip(names, all_boxes); print(" $(n)=$(t)"); end; println()

    for (name, tag) in zip(names, all_boxes)
        faces = gmsh.model.getBoundary([(3,tag)], false, false, false)
        nf = length(faces)
        println("  $(name) (vol $(tag)): $(nf) faces $(nf==6 ? "✓" : "✗ WARNING")")
        nf == 6 || error("Volume $(name) has $(nf) faces, expected 6")
    end

    # ── Node counts ────────────────────────────────────────────────────────────
    # dy_hi controls resolution in tall fluid region — keep coarse (2-4mm)
    nx_left   = max(2, round(Int, (x1-x0) / fin_pitch) + 1)
    nx_porous = max(2, round(Int, (x2-x1) / fin_pitch) + 1)
    nx_right  = max(2, round(Int, (x3-x2) / fin_pitch) + 1)
    ny_lo     = max(2, round(Int, (y1-y0) / fin_pitch) + 1)
    ny_hi     = max(2, round(Int, (y2-y1) / dy_hi)     + 1)
    nz_nodes  = nz + 1

    println("Nodes — X: left=$(nx_left) porous=$(nx_porous) right=$(nx_right)")
    println("        Y: lo=$(ny_lo) hi=$(ny_hi)  Z=$(nz_nodes)")
    println("Estimated hex count: $(nx_porous * ny_lo * nz) (porous) + approx fluid")

    tol = 1e-6

    for (_, ctag) in gmsh.model.getEntities(1)
        xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(1, ctag)
        dx = abs(xmax-xmin); dy = abs(ymax-ymin); dz = abs(zmax-zmin)

        n = if dx > tol && dy < tol && dz < tol
            span = (round(min(xmin,xmax),digits=4), round(max(xmin,xmax),digits=4))
            if     abs(span[1]-x0)<tol && abs(span[2]-x1)<tol; nx_left
            elseif abs(span[1]-x1)<tol && abs(span[2]-x2)<tol; nx_porous
            elseif abs(span[1]-x2)<tol && abs(span[2]-x3)<tol; nx_right
            else; 3
            end
        elseif dy > tol && dx < tol && dz < tol
            span = (round(min(ymin,ymax),digits=4), round(max(ymin,ymax),digits=4))
            if     abs(span[1]-y0)<tol && abs(span[2]-y1)<tol; ny_lo
            elseif abs(span[1]-y1)<tol && abs(span[2]-y2)<tol; ny_hi
            else; 3
            end
        else
            nz_nodes
        end

        gmsh.model.mesh.setTransfiniteCurve(ctag, n)
    end

    for tag in all_boxes
        for (_, s) in gmsh.model.getBoundary([(3,tag)], false, false, false)
            gmsh.model.mesh.setTransfiniteSurface(abs(s))
            gmsh.model.mesh.setRecombine(2, abs(s))
        end
        gmsh.model.mesh.setTransfiniteVolume(tag)
    end

    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    # ── Physical groups ────────────────────────────────────────────────────────
    fluid_vols = [b_left_lo, b_left_hi, b_above, b_right_lo, b_right_hi]
    gmsh.model.addPhysicalGroup(3, fluid_vols, -1, "fluid")
    gmsh.model.addPhysicalGroup(3, [b_porous],  -1, "fins")

    inlet_s = Int[]; outlet_s = Int[]
    for (dim, tag) in gmsh.model.getEntities(2)
        _,_,zmin,_,_,zmax = gmsh.model.getBoundingBox(dim, tag)
        if   zmin > z1-0.01 && zmax > z1-0.01; push!(inlet_s,  tag)
        elseif zmin < 0.01  && zmax < 0.01;    push!(outlet_s, tag)
        end
    end
    gmsh.model.addPhysicalGroup(2, inlet_s,  -1, "inlet")
    gmsh.model.addPhysicalGroup(2, outlet_s, -1, "outlet")

    fluid_surfs  = Set(abs(tag) for (_,tag) in
        Iterators.flatten(gmsh.model.getBoundary([(3,v)], false, false, false) for v in fluid_vols))
    porous_surfs = Set(abs(tag) for (_,tag) in
        gmsh.model.getBoundary([(3,b_porous)], false, false, false))
    interface    = collect(intersect(fluid_surfs, porous_surfs))
    gmsh.model.addPhysicalGroup(2, interface, -1, "fluid_porous_interface")

    all_surfs  = Set(tag for (_,tag) in gmsh.model.getEntities(2))
    used_surfs = Set(vcat(inlet_s, outlet_s, interface))
    gmsh.model.addPhysicalGroup(2, collect(setdiff(all_surfs, used_surfs)), -1, "walls")

    gmsh.model.mesh.generate(3)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    hex_tags, _     = gmsh.model.mesh.getElementsByType(5)
    tet_tags, _     = gmsh.model.mesh.getElementsByType(4)
    all_elems       = vcat(hex_tags, tet_tags)

    println("\n── Mesh Statistics ──────────────────────")
    println("  Nodes: $(length(node_tags))")
    println("  Hexes: $(length(hex_tags))  Tets: $(length(tet_tags))")
    if !isempty(all_elems)
        q = gmsh.model.mesh.getElementQualities(all_elems, "minSICN")
        println("  Min SICN:  $(round(minimum(q), digits=4))")
        println("  Mean SICN: $(round(sum(q)/length(q), digits=4))")
    end
    println("─────────────────────────────────────────")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(outfile)
    gmsh.finalize()
    println("Written: $outfile")
end

build_mesh(4.74, 1.60, 0.23, outfile="Julia_Attempt/coldplate.msh", nz=10, dy_hi=3.0)