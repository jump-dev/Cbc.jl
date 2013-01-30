@unix_only begin
if OS_NAME == :Linux
    shlib_ext = ".so"
elseif OS_NAME == :Darwin
    shlib_ext = ""
else
    error("Platform not currently supported")
end
const coinmp_lib = joinpath(Pkg.dir(),"CoinMP","deps","usr","lib","libCoinMP$shlib_ext")
end

@windows_only begin
const coinmp_lib = joinpath(Pkg.dir(),"CoinMP","deps","CoinMP.dll")
end

