using BinDeps

@BinDeps.setup

@unix_only begin
    libclp = library_dependency("libclp",aliases=["libClp"])
    libcoinmp = library_dependency("libcoinmp",aliases=["libCoinMP"])
end
@windows_only begin
    if Int != Int32
        error("Win64 platform is not yet supported by this package")
    end
    libclp = library_dependency("libclp",aliases=["CoinMP"])
    libcoinmp = library_dependency("libcoinmp",aliases=["CoinMP"])
end

coinmpname = "CoinMP-1.7.6"

provides(Sources, URI("http://www.coin-or.org/download/source/CoinMP/$coinmpname.tgz"),
    [libclp,libcoinmp], os = :Unix)

provides(Binaries, URI("http://www.mit.edu/~mlubin/CoinMP_julia_20130903.tar.gz"),
    [libclp,libcoinmp], os = :Windows)

@osx_only begin
    using Homebrew
    provides( Homebrew.HB, "coinmp", [libclp, libcoinmp], os = :Darwin )
end

prefix=joinpath(BinDeps.depsdir(libclp),"usr")
patchdir=BinDeps.depsdir(libclp)
srcdir = joinpath(BinDeps.depsdir(libclp),"src",coinmpname) 

provides(SimpleBuild,
    (@build_steps begin
        GetSources(libclp)
        @build_steps begin
            ChangeDirectory(srcdir)
            `cat $patchdir/CoinMP-emptyproblem.patch` |> `patch -N -p1`
            `./configure --prefix=$prefix --enable-dependency-linking`
            `make install`
        end
    end),[libclp,libcoinmp], os = :Unix)

@BinDeps.install [:libclp => :libclp, :libcoinmp => :libcoinmp]
