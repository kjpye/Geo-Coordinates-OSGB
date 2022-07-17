unit class Geo::Coordinates::OSGB;

use Geo::Geometry;

my constant \deg2rad = π ÷ 180;

my \N0 = -100000;   # northing of true origin
my \E0 = 400000;   # easting of true origin
my \F0 = 0.9996012717;   # scale factor on central meridian
my \φ0 = 49 × deg2rad;       # latitude of true origin
my \λ0 = -2.0 × deg2rad;     # longitude of true origin

# ellipsoid
my \a = 6377563.396; # Airy
my \b = 6356256.909;
my \e2 = (a² - b²) ÷ a²;
my \n = (a - b) ÷ (a + b);
my \f = (a - b) ÷ a;

multi sub latlon-to-osgb (Point $latlon) {
    samewith($latlon.y, $latlon.x);
}

multi sub latlon-to-osgb (:$latitude, :$longitude) {
    samewith($latitude, $longitude);
}

multi sub latlon-to-osgb ($latitude, $longitude) {
    my \φ = $latitude × deg2rad;
    my \λ = $longitude × deg2rad;
    my \ν = a × F0 ÷ sqrt(1 - e2 × sin(φ)²);
    my \ρ = a × F0 × (1 - e2) ÷ (1 - e2 × sin(φ)²) ** 1.5;
    my \η2 = ν ÷ ρ - 1;
    my \M = b × F0 × (  (1 + n + 5/4 × n² + 5/4 × n³) × (φ - φ0)
                      - (3×n + 3×n² + 21/8 × n³) × sin(φ - φ0) × cos(φ + φ0)
                      + (15/8 × n² + 15/8 × n³) × sin(2×(φ - φ0)) × cos(2×(φ + φ0))
                      - (35/24 × n³ × sin(3×(φ - φ0)) × cos(3×(φ + φ0)))
                     );
    my \I    = M + N0;
    my \II   = ν/2 × sin(φ) × cos(φ);
    my \III  = ν/24 × sin(φ) × cos(φ)³×(5 - tan(φ)² + 9 × η2);
    my \IIIA = ν/720 × sin(φ) × cos(φ)⁵×(61 - 58×tan(φ)² + tan(φ)⁴);
    my \IV   = ν × cos(φ);
    my \V    = ν/6 × cos(φ)³ × (ν ÷ ρ - tan(φ)²);
    my \VI   = ν/120 × cos(φ)⁵ × (5 - 18 × tan(φ)² + tan(φ)⁴ + 14 × η2 - 58 × tan(φ)² × η2);
    my \N = I + II × (λ - λ0)² + III × (λ - λ0)⁴ + IIIA × (λ - λ0)⁶;
    my \E = E0 + IV × (λ - λ0) + V × (λ - λ0)³ + VI × (λ - λ0)⁵;
    Point.new(x => E, y => N);
}

dd latlon-to-osgb(52.6575703056×deg2rad, 1.71792158333×deg2rad);
dd latlon-to-osgb(Point.new(y => 52.6575703056×deg2rad, x => 1.71792158333×deg2rad));

multi sub osgb-to-latlon (Point $point) {
    samewith($point.y, $point.x);
}

multi sub osgb-to-latlon(:$easting, :$northing) {
    samewith($easting, $northing);
}

multi sub osgb-to-latlon (Real \E, Real \N) {
    my $M = 0;
    my $φ = φ0;
    while N - N0 - $M ≥ 0.00001 {
         $φ = (N - N0 - $M) ÷ (a × F0) + $φ;
      my $ν = a × F0 ÷ sqrt(1 - e2 × sin($φ)²);
      my $ρ = a × F0 × (1 - e2) ÷ (1 - e2 × sin($φ)²) ** 1.5;
      my $η2 = $ν ÷ $ρ - 1;
         $M = b × F0 × (  (1 + n + 5/4 × n² + 5/4 × n³) × ($φ - φ0)
                        - (3×n + 3×n² + 21/8 × n³) × sin($φ - φ0) × cos($φ + φ0)
                        + (15/8 × n² + 15/8 × n³) × sin(2×($φ - φ0)) × cos(2×($φ + φ0))
                        - (35/24 × n³ × sin(3×($φ - φ0)) × cos(3×($φ + φ0)))
                       );
    }
    my \ν = a × F0 ÷ sqrt(1 - e2 × sin($φ)²);
    my \ρ = a × F0 × (1 - e2) ÷ (1 - e2 × sin($φ)²) ** 1.5;
    my \η2 = ν ÷ ρ - 1;
    my \VII = tan($φ) ÷ (2 × ρ × ν);
    my \VIII = tan($φ) ÷ (24 × ρ × ν³) × (5 + 3 × tan($φ)² + η2 - 9 × tan($φ)² × η2);
    my \IX = tan($φ) ÷ (720 × ρ × ν⁵) × (61 + 90 × tan($φ)² + 45 × tan($φ)⁴);
    my \X = sec($φ) ÷ ν;
    my \XI = sec($φ) ÷ (6 × ν³) × (ν/ρ + 2 × tan($φ)²);
    my \XII = sec($φ) ÷ (120 × ν⁵) × (5 + 28 × tan($φ)² + 24 × tan($φ)⁴);
    my \XIIA = sec($φ) ÷ (5040 × ν⁷) × (61 + 662 × tan($φ)² + 1320 × tan($φ)⁴ + 720 × tan($φ)⁶);
    my \φ = ($φ - VII × (E - E0)² + VIII × (E - E0)⁴ - IX × (E - E0)⁶) ÷ deg2rad;
    my \λ = (λ0 + X × (E - E0) - XI × (E - E0)³ + XII × (E -E0)⁵ - XIIA × (E - E0)⁷) ÷ deg2rad;
    Point.new(x => λ, y => φ);
}

dd osgb-to-latlon(651409.9029626573, 313177.2703357119);
dd osgb-to-latlon(Point.new(313177.2703357119, 651409.9029626573));

=begin pod

=head1 NAME

Geo::Coordinates::OSGB - Convert coordinates between the Ordnance Survey National Grid and latitude and longitude

=head1 SYNOPSIS

=begin code :lang<raku>

use Geo::Geometry;
use Geo::Coordinates::OSGB;

my $point1 = latlon-to-usgb(1.72, 52.65);
my $point2 = usgb-to-latlon($point1);

=end code

=head1 DESCRIPTION

Geo::Coordinates::OSGB provides routines to convert coordinates between the Ordnance Survey National Grid and latitude and longitude.

The routines used are based on Appendix C of I<A Guide to Coordinate Systems in Great Britain>, published by the Ordnance Survey.

All latitudes and longitudes are expressed in degrees, and eastings and northings in metres.

Every routine returns a Geo::Geometry::Point object.

Two routines are available, each of which takes its arguments in various forms:

=head2 latlon-to-usgb

C<latlon-to-usgb> takes a point expressed as a latitude and a longitude (based on the Airy ellipsoid) and converts it to a C<Geo::Geometry::Point> object which contains the easting and northing of that point on the Ordnance Survey grid.

The arguments can be given as a C<Point> object, as named arguments C<latitude> and C<longitude>, or as unnamed arguments, with the latitude given first. The following three calls to C<latlon-to-usgb> are equivalent:

=begin code :lang<raku>
    latlon-to-usgb(52, 1.3);
    latlon-to-usgb(longitude => 1.3, latitude => 52);
    latlon-to-usgb(Point.new(x => 1.3. y => 52));
=end code
    
=head2 usgb-to-latlon
=head1 AUTHOR

Kevin Pye <kjpraku@pye.id.au>

=head1 COPYRIGHT AND LICENSE

Copyright 2022 Kevin Pye

This library is free software; you can redistribute it and/or modify it under the Artistic License 2.0.

=end pod
