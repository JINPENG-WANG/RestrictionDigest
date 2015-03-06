use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'RestrictionDigest',
    AUTHOR           => q{Jinpeng Wang, Li Li, Haigang Qi, Xuedi Du, and Guofan Zhang <restrictiondigest@163.com>},
    VERSION_FROM     => 'lib/RestrictionDigest.pm',
    ABSTRACT_FROM    => 'lib/RestrictionDigest.pm',
    LICENSE          => 'Artistic_2_0',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        'IO::File'              => 1.14,
        'File::Basename' => 2.74,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'RestrictionDigest-*' },
);
