include("CHT_Functions.jl")

stp_dir = "/home/aaronbest/Desktop/GodelaChallenge/coldplate"
prefix = "coldplate_inst"
suffix = ".stp"

N_runs = 3

for i in 1:N_runs

    #Setup directory and locate files
    if !isdir("Case"* "_" *string(i-1))
        mkdir("Case"* "_" *string(i-1))
    end
    case_dir = "/home/aaronbest/Desktop/GodelaChallenge/Godela_Code"* "/" *"Case"* "_" *string(i-1)
    input_stp = stp_dir* "/" *prefix*string(i-1)*suffix
    
    #Create Mesh
    if !isfile(case_dir* "/" *"coldplate.msh")
        print("Creating Mesh...")
        CreateMesh(case_dir,input_stp,false)
        print("Mesh Done...")
    else
        print("Mesh file found, skipping step...")
    end

    #Solve CHT
    print("Solving CHT...")
    SolveCHT(case_dir)
    print("CHT Solved...")
end
