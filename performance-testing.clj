(def n 2000) ; Matrix size
(def k 1000) ; Index of the returned element - force lazy sequence to be evaluated.
(def A (take n (repeatedly (fn [] (take n (repeatedly #(rand 10)))))))
;(def v (take n (repeatedly #(rand 10))))
(def v (range n))

(time (nth (map #(apply + (map * % v)) A) k))
(time (nth (map (fn [ai] (reduce #(+ %1 (apply * %2)) 0 (map vector v ai))) A) k))
