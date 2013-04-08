using BinDeps

@unix_only begin
libname = "CoinMP-1.6.0"
prefix = joinpath(Pkg.dir(),"CoinMP","deps","usr")
libdir = joinpath(JULIA_HOME,"..","lib")
if !isfile("$libname.tgz")
    run(download_cmd("http://www.coin-or.org/download/source/CoinMP/$libname.tgz","$libname.tgz"))
    run(`tar xvzf $libname.tgz`)
    cd("$libname")
    run(`cat ../CoinMP-makefile.patch` | `patch -p1`)
    run(`cat ../CoinMP-strcmp.patch` | `patch -p1`)
    run(`./configure --prefix=$prefix`)
    run(`make install`)
end
end # unix_only

@windows_only begin
    if !isfile("CoinMP_julia.tar.gz")
        run(download_cmd("http://www.mit.edu/~mlubin/CoinMP_julia.tar.gz","CoinMP_julia.tar.gz"))
    end
    if !isfile("CoinMP.dll")
        run(unpack_cmd("CoinMP_julia.tar.gz","."))
    end
end


