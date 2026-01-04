package api

import (
	"net/http"
	"net/http/httptest"
	"net/url"
	"reflect"
	"strings"
	"testing"
)

func TestParseCategoryFilter(t *testing.T) {
	t.Run("absent", func(t *testing.T) {
		filter, ok := parseCategoryFilter(url.Values{})
		if ok {
			t.Fatalf("expected ok=false, got true")
		}
		if filter != nil {
			t.Fatalf("expected nil filter, got %#v", filter)
		}
	})

	t.Run("commaSeparated", func(t *testing.T) {
		q, _ := url.ParseQuery("categories=T,B")
		filter, ok := parseCategoryFilter(q)
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})

	t.Run("jsonArray", func(t *testing.T) {
		q, _ := url.ParseQuery(`categories=["T","B"]`)
		filter, ok := parseCategoryFilter(q)
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})

	t.Run("jsonEmpty", func(t *testing.T) {
		q, _ := url.ParseQuery(`categories=[]`)
		filter, ok := parseCategoryFilter(q)
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		if filter == nil || len(filter) != 0 {
			t.Fatalf("expected non-nil empty filter, got %#v", filter)
		}
	})

	t.Run("emptyString", func(t *testing.T) {
		q, _ := url.ParseQuery(`categories=`)
		filter, ok := parseCategoryFilter(q)
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		if filter == nil || len(filter) != 0 {
			t.Fatalf("expected non-nil empty filter, got %#v", filter)
		}
	})

	t.Run("repeatedParams", func(t *testing.T) {
		q := url.Values{"categories": {"T", "B"}}
		filter, ok := parseCategoryFilter(q)
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})
}

func TestParseCategoryFilterBody(t *testing.T) {
	t.Run("empty", func(t *testing.T) {
		r := httptest.NewRequest(http.MethodPost, "/d/default/tiles/0/0/0/category/cell_type.png", strings.NewReader(""))
		filter, ok, err := parseCategoryFilterBody(r)
		if err != nil {
			t.Fatalf("expected err=nil, got %v", err)
		}
		if ok {
			t.Fatalf("expected ok=false, got true")
		}
		if filter != nil {
			t.Fatalf("expected nil filter, got %#v", filter)
		}
	})

	t.Run("jsonArray", func(t *testing.T) {
		r := httptest.NewRequest(http.MethodPost, "/d/default/tiles/0/0/0/category/cell_type.png", strings.NewReader(`["T","B"]`))
		filter, ok, err := parseCategoryFilterBody(r)
		if err != nil {
			t.Fatalf("expected err=nil, got %v", err)
		}
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})

	t.Run("jsonObject", func(t *testing.T) {
		r := httptest.NewRequest(http.MethodPost, "/d/default/tiles/0/0/0/category/cell_type.png", strings.NewReader(`{"categories":["T","B"]}`))
		filter, ok, err := parseCategoryFilterBody(r)
		if err != nil {
			t.Fatalf("expected err=nil, got %v", err)
		}
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})

	t.Run("formEncodedJson", func(t *testing.T) {
		r := httptest.NewRequest(http.MethodPost, "/d/default/tiles/0/0/0/category/cell_type.png", strings.NewReader(`categories=["T","B"]`))
		filter, ok, err := parseCategoryFilterBody(r)
		if err != nil {
			t.Fatalf("expected err=nil, got %v", err)
		}
		if !ok {
			t.Fatalf("expected ok=true, got false")
		}
		want := []string{"T", "B"}
		if !reflect.DeepEqual(filter, want) {
			t.Fatalf("expected %#v, got %#v", want, filter)
		}
	})
}
