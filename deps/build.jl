using BinDeps

@BinDeps.setup

@unix_only begin
    libclp = library_dependency("libclp",aliases=["libClp"])
    libcbcsolver = library_dependency("libcbcsolver",aliases=["libCbcSolver"])
end
@windows_only begin
    using WinRPM
    push!(WinRPM.sources, "http://download.opensuse.org/repositories/home:/kelman:/mingw-coinor/openSUSE_13.1")
    libclp = library_dependency("libclp",aliases=["libClp-1"])
    libcbcsolver = library_dependency("libcbcsolver",aliases=["libCbcSolver-3"])
    provides(WinRPM.RPM, "coin-or-Cbc", [libclp,libcbcsolver], os = :Windows)
end

coinmpname = "CoinMP-1.7.6"

provides(Sources, URI("http://www.coin-or.org/download/source/CoinMP/$coinmpname.tgz"),
    [libclp,libcbcsolver], os = :Unix)

@osx_only begin
    using Homebrew
    provides( Homebrew.HB, "coinmp", [libclp, libcbcsolver], os = :Darwin )
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
    end),[libclp,libcbcsolver], os = :Unix)

@BinDeps.install [:libclp => :libclp, :libcbcsolver => :libcbcsolver]
