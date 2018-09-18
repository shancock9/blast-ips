package Blast::IPS::Utils;

# MIT License
# Copyright (c) 2018 Steven Hancock
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Some common utility function

use strict;
use warnings;
use 5.006;

my $VERSION = 1.00;
use Carp;

our @EXPORT_OK = qw(
  check_keys
);
use Exporter;
our @ISA = qw(Exporter);

sub check_keys {
    my ( $rtest, $rvalid, $msg, $exact_match ) = @_;

    # Check the keys of a hash for validity:
    # $rtest  = ref to hash to test
    # $rvalid = ref to has with valid keys

    # $msg = a message to write in case of error
    # $exact_match defines the type of check:
    #     = false: test hash must not have unknown key
    #     = true:  test hash must have exactly same keys as known hash
    my @unknown_keys =
      grep { $_ && !exists $rvalid->{$_} } keys %{$rtest};
    my @missing_keys =
      grep { !exists $rtest->{$_} } keys %{$rvalid};
    my $error = @unknown_keys;
    if ($exact_match) { $error ||= @missing_keys }
    if ($error) {
        local $" = ')(';
        my @expected_keys = sort keys %{$rvalid};
        @missing_keys = sort @missing_keys;
        @unknown_keys = sort @unknown_keys;
	my $caller = caller();
        croak(<<EOM);
------------------------------------------------------------------------
$caller program error detected checking hash keys
Message is: '$msg'
Valid keys are: (@expected_keys)
Keys not seen : (@missing_keys)
Unknown key(s): (@unknown_keys)
------------------------------------------------------------------------
EOM
    }
    return;
}

1;
