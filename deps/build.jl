using BinDeps

@BinDeps.setup

@unix_only begin
    libclp = library_dependency("libclp",aliases=["libClp"])
    libcoinmp = library_dependency("libcoinmp",aliases=["libCoinMP"])
end
@windows_only begin
    libclp = library_dependency("libclp",aliases=["CoinMP"])
    libcoinmp = library_dependency("libcoinmp",aliases=["CoinMP"])
end

coinmpname = "CoinMP-1.7.0"

provides(Sources, URI("http://www.coin-or.org/download/source/CoinMP/$coinmpname.tgz"),
    [libclp,libcoinmp], os = :Unix)

# TODO: Fix this tarball
provides(Binaries, URI("http://www.mit.edu/~mlubin/CoinMP_julia.tar.gz"),
    [libclp,libcoinmp], os = :Windows)

prefix=joinpath(BinDeps.depsdir(libclp),"usr")
patchdir=BinDeps.depsdir(libclp)
srcdir = joinpath(BinDeps.depsdir(libclp),"src",coinmpname) 

provides(SimpleBuild,
    (@build_steps begin
        GetSources(libclp)
        @build_steps begin
            ChangeDirectory(srcdir)
            `cat $patchdir/CoinMP-makefile.patch` |> `patch -p1`
            `cat $patchdir/CoinMP-strcmp.patch` |> `patch -p1`
            `cat $patchdir/CoinMP-loglevel.patch` |> `patch -p1`
            `cat $patchdir/Clp-interface.patch` |> `patch -p0`
            `./configure --prefix=$prefix`
            `make install`
        end
    end),[libclp,libcoinmp], os = :Unix)

@BinDeps.install
