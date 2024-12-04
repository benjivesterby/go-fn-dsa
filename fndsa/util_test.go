package fndsa

import (
	"testing"
)

func TestSigningKeySize(t *testing.T) {
	var expected = [11]int{
		0, 0, 13, 25, 49, 97, 177, 353, 641, 1281, 2305,
	}
	for logn := uint(2); logn <= 10; logn++ {
		s := SigningKeySize(logn)
		if s != expected[logn] {
			t.Fatalf("ERR: logn=%d -> %d (exp: %d)\n", logn, s, expected[logn])
		}
	}
}

func TestVerifyingKeySize(t *testing.T) {
	var expected = [11]int{
		0, 0, 8, 15, 29, 57, 113, 225, 449, 897, 1793,
	}
	for logn := uint(2); logn <= 10; logn++ {
		s := VerifyingKeySize(logn)
		if s != expected[logn] {
			t.Fatalf("ERR: logn=%d -> %d (exp: %d)\n", logn, s, expected[logn])
		}
	}
}

func TestSignatureSize(t *testing.T) {
	var expected = [11]int{
		0, 0, 47, 52, 63, 82, 122, 200, 356, 666, 1280,
	}
	for logn := uint(2); logn <= 10; logn++ {
		s := SignatureSize(logn)
		if s != expected[logn] {
			t.Fatalf("ERR: logn=%d -> %d (exp: %d)\n", logn, s, expected[logn])
		}
	}
}
