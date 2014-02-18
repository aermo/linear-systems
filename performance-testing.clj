(load-file "numerical-methods.clj")
(ns
  ^{:author "Antti Rämö",
     :doc ""}
  numerical-methods-performance
  (:use [numerical-methods :as nm]))

(def n 2000) ; Matrix size
(def k 1000) ; Index of the returned element - force lazy sequence to be evaluated.
(def A (take n (repeatedly (fn [] (take n (repeatedly #(rand 10)))))))
(def v (take n (repeatedly #(rand 10))))
;(def v (range n))

;; Random sequences are evaluated here.
;(time (nth (map #(apply + (map * % v)) A) k))

;; Performance comparison.
;(time (nth 
;  (map #(apply + (map * % v)) A)
;  k))
;(time (nth 
;  (map (fn [ai] (reduce #(+ %1 (apply * %2)) 0 (map vector v ai))) A)
;  k))
;; --> Simple is better!




;(time (nth (nm/backward-substitution A v) k))

;(time (nth (nm/backward-substitution A v) k))
;(time (nth (nm/backward-substitution1 A v) k))
;(time (nth (nm/backward-substitution2 A v) k))

(def n 300) ; Matrix size
(def A (take n (repeatedly (fn [] (take n (repeatedly #(rand 10)))))))
(def v (take n (repeatedly #(rand 10))))

(def res1a (time (nm/lu-factorization-sparse A)))
(def res1b (time (nm/lu-factorization-sparse A)))
(def res2a (time (nm/lu-factorization-sparse2 A)))
(def res2b (time (nm/lu-factorization-sparse2 A)))
(def res3a (time (nm/lu-factorization-sparse3 A)))
(def res3b (time (nm/lu-factorization-sparse3 A)))
(def res1a2 (time (nm/lu-factorization-sparse A)))
(def res1b2 (time (nm/lu-factorization-sparse A)))
(def res2a2 (time (nm/lu-factorization-sparse2 A)))
(def res2b2 (time (nm/lu-factorization-sparse2 A)))
(def res3a2 (time (nm/lu-factorization-sparse3 A)))
(def res3b2 (time (nm/lu-factorization-sparse3 A)))
